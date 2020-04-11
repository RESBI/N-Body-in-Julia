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
# ⁽ᶜᵒᵈᵉ ᵇʸ ᴿᵉˢᵇⁱ⁾

# Modules
using LinearAlgebra

# Constants
const n = 1024*16
const G = 6.67408e-11
const Δt = 60*12
const Δt² = Δt^2
const Δtₕ = Δt/2
const c = 299792458 # m/s
const Δt⁻²c⁻² = (c*Δt)^(-2)
const Mₘₐₓ = 1e22
const Mₘᵢₙ = 1e20
const Areaₘᵢₙ = [-1e11 -1e11 -1e11]
const Areaₘₐₓ = [1e11 1e11 1e11]
const Areaₜₑₘₚ = Areaₘₐₓ - Areaₘᵢₙ
const v⃗ₘᵢₙ  = [-5e2 -5e2 -5e2]
const v⃗ₘₐₓ  = [5e2 5e2 5e2]
const v⃗ₜₑₘₚ  = v⃗ₘₐₓ  - v⃗ₘᵢₙ 
const SavePer = 5
const Total = 1024*16

# Generate
MGΔtΔtₕ = (rand(n)*(Mₘₐₓ-Mₘᵢₙ) + ones(n)*Mₘᵢₙ)*G*Δt*Δtₕ
P₀ = zeros(3, n)
P₁ = zeros(3, n)
Δtv⃗ₕ  = zeros(3, n)
#Δt²a⃗ₕ = zeros(3, n)
for i=1:3
  P₀[i, :] = rand(n)*Areaₜₑₘₚ[i] + ones(n)*Areaₘᵢₙ[i]
  Δtv⃗ₕ[i,  :] = (rand(n)*v⃗ₜₑₘₚ[i]  + ones(n)*v⃗ₘᵢₙ[i])*Δtₕ
end

#Set a blackhole
P₀[:, 1] = [0.0 0.0 0.0]
Δtv⃗ₕ[:, 1] = [0.0 0.0 0.0]
const MGΔtΔtₕ[1] = 1e30*G*Δt*Δtₕ

# Functions
function update(SYNC_ARRAY, P₀, P₁, Δtv⃗ₕ)
  @inbounds Threads.@threads for i=1:n
  #for i=1:n
    Δt²a⃗ₕ  = [0.0, 0.0, 0.0]
    #γ = sqrt(1 - (norm(Δtv⃗ₕ[i])^2)*Δt⁻²c⁻²)
    for j=1:(i-1)
      Δt²a⃗ₕ  += MGΔtΔtₕ[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    for j=(i+1):n
      Δt²a⃗ₕ  += MGΔtΔtₕ[j]/(norm(P₀[:, j] - P₀[:, i])^3)*(P₀[:, j] - P₀[:, i])
    end
    P₁[:, i] = P₀[:, i] + 2*Δtv⃗ₕ[:, i] + Δt²a⃗ₕ #*γ
    Δtv⃗ₕ[:, i] = Δtv⃗ₕ[:, i] + Δt²a⃗ₕ #*γ
    # Mark
    #SYNC_ARRAY[i] = 1
  end
end

# Run
for j=0:Total-1
  SYNC_ARRAY = zeros(n)
  if mod(j, 2)==0
    @time update(SYNC_ARRAY, P₀, P₁, Δtv⃗ₕ)
  else
    @time update(SYNC_ARRAY, P₁, P₀, Δtv⃗ₕ)
  end
  # Wait
  #while sum(SYNC_ARRAY) < n
  #  pass
  #end
  if mod(j, SavePer)==0
    println(Δtv⃗ₕ[:, 1]/Δtₕ)
    println(P₀[:, 1])
  end
end
