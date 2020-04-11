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
Δt² = Δt^2
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
MGΔtΔtₕ = (rand(n)*(Mₘₐₓ-Mₘᵢₙ) + ones(n)*Mₘᵢₙ)*G*Δt*Δtₕ
P₀ = zeros(3, n)
P₁ = zeros(3, n)
Δtv⃗ₕ  = zeros(3, n)
for i=1:3
  P₀[i, :] = rand(n)*Areaₜₑₘₚ[i] + ones(n)*Areaₘᵢₙ[i]
  Δtv⃗ₕ[i,  :] = (rand(n)*v⃗ₜₑₘₚ[i]  + ones(n)*v⃗ₘᵢₙ[i])*Δtₕ
end

#Set an blackhole
P₀[:, 1] = [0.0 0.0 0.0]
Δtv⃗ₕ[:, 1] = [0.0 0.0 0.0]
MGΔtΔtₕ[1] = 1e30*G*Δt*Δtₕ

# Functions
#rsqrt(x) = 1/sqrt(x);
function update(P₀, P₁, Δtv⃗ₕ)
  for i=1:n
    Δt²a⃗ₕ  = [0.0, 0.0, 0.0]
    for j=1:(i-1)
      Δt²a⃗ₕ  += MGΔtΔtₕ[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    for j=(i+1):n
      Δt²a⃗ₕ  += MGΔtΔtₕ[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    P₁[:, i] = P₀[:, i] + 2*Δtv⃗ₕ[:, i] + Δt²a⃗ₕ
    Δtv⃗ₕ[:, i] = Δtv⃗ₕ[:, i] + Δt²a⃗ₕ
  end
end

# Run
for j=0:Total-1
  @time if mod(j, 2)==0
    update(P₀, P₁, Δtv⃗ₕ)
  else
    update(P₁, P₀, Δtv⃗ₕ)
  end
  if mod(j, SavePer)==0
    println(Δtv⃗ₕ[:, 1]/Δtₕ)
    println(P₀[:, 1])
  end
end
