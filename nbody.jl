# We have n particles here.
# For each particle i (i ∈ N⁺)  :
# We have Mₘᵢₙ ≤ mᵢ ≤ Mₘₐₓ
# We have an three dimentional vector v⃗ᵢ for i's velocity. 
# And an three dimentional point Pᵢ for i's locationn.
# We have Gravitational Constant G ≈ 6.67408*10⁻¹¹ (m³kg⁻¹s⁻²) .
# And a small time step Δt→0 (actually, it's a very small value).
# At the beginning, every particle will spread in an area randomly.
# At the beginning, every particle will have a velocity which direction and size are all given randomly but in a range.
# 

# Modules
using LinearAlgebra

# Constants
n = 1000
G = 6.67408e-11
Δt = 60*12
Δtₕ = Δt/2
Mₘₐₓ = 1e22
Mₘᵢₙ = 1e20
Areaₘᵢₙ = [-1e11 -1e11 -1e11]
Areaₘₐₓ = [1e11 1e11 1e11]
Areaₜₑₘₚ = Areaₘₐₓ - Areaₘᵢₙ
v⃗ₘᵢₙ  = [-5e2 -5e2 -5e2]
v⃗ₘₐₓ  = [5e2 5e2 5e2]
v⃗ₜₑₘₚ  = v⃗ₘₐₓ  - v⃗ₘᵢₙ 
SavePer = 5
Total = 1024*16

# Generate
MGΔt = (rand(n)*(Mₘₐₓ-Mₘᵢₙ) + ones(n)*Mₘᵢₙ)*G*Δt
P₀ = zeros(3, n)
P₁ = zeros(3, n)
v⃗  = zeros(3, n)
for i=1:3
  P₀[i, :] = rand(n)*Areaₜₑₘₚ[i] + ones(n)*Areaₘᵢₙ[i]
  v⃗[i,  :] = rand(n)*v⃗ₜₑₘₚ[i]  + ones(n)*v⃗ₘᵢₙ[i]
end

#Set an blackhole
P₀[:, 1] = [0.0 0.0 0.0]
v⃗[:, 1] = [0.0 0.0 0.0]
MGΔt[1] = 1e30*G*Δt

# Functions
#rsqrt(x) = 1/sqrt(x);
function update(P₀, P₁, v⃗)
  for i=1:n
    Δta⃗  = [0.0, 0.0, 0.0]
    for j=1:(i-1)
      Δta⃗  += MGΔt[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    for j=(i+1):n
      Δta⃗  += MGΔt[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    P₁[:, i] = P₀[:, i] + v⃗[:, i]*Δt + Δta⃗ *Δtₕ
    v⃗[:, i] = v⃗[:, i] + Δta⃗
  end
end

# Run
for j=0:Total-1
  @time if mod(j, 2)==0
    update(P₀, P₁, v⃗)
  else
    update(P₁, P₀, v⃗)
  end
  if mod(j, SavePer)==0
    println(P₀[:, 1])
  end
end
