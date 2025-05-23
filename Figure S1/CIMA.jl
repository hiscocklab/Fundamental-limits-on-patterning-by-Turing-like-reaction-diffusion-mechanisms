## import and load
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics
include("required_scripts_2.jl")

## Specify the model
model = @reaction_network begin
    k₁,             ∅ --> A
    1.0,            A --> ∅
    4.0*I/(1+A^2),  A --> ∅
    k₂*A,           ∅ --> I
    k₂*A/(1+A^2),   I --> ∅
end

## Specify model parameters
params = model_parameters()
Nscreen = 100

params.reaction["k₁"] = screen_values(min = 10, max = 30, mode = "linear", number = N_screen)
params.reaction["k₂"] = screen_values(min = 0.1, max = 30, mode = "linear", number = N_screen)

params.diffusion["A"] = [1.0] 
params.diffusion["I"] = [15.0]


## Compile functions and solvers

include("simulate_scripts_2.jl")

## Testing for Turing

@time turing_params = returnTuringParams(model, params);

## Plotting 

using Plots
U₀ = VectorOfArray(turing_params.steady_state_values)[1,:]
V₀ = VectorOfArray(turing_params.steady_state_values)[2,:]
a = get_param(model, turing_params,"k₁","reaction")
b = get_param(model, turing_params,"k₂","reaction")
Dᵥ = get_param(model, turing_params,"I","diffusion")

ã = range(10,stop=30,length=30)
Ũ₀ = ã / 5
Ṽ₀ = 1 .+ã.^2 ./ 25.
D̃ᵥ = Dᵥ[1]
i = findall(Dᵥ .== D̃ᵥ)
b̃H = 0.6*ã .- 25 ./ ã
b̃T = (D̃ᵥ ./ (5*ã)) .* (13*ã.^2 .+ 125 .- 4*ã .* sqrt.(10*(ã.^2 .+ 25)))

Plots.plot(ã,[b̃H,b̃T], labels=["lower bound" "upper bound"],xlabel="b", ylabel="a",linewidth = 5)

Plots.scatter!(a[i],b[i],strokewidth = 0, label = "numerical sol.", color = :black, markersize=1)
title!("Analytical and numerical solution for CIMA model")
