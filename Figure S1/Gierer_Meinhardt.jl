## import and load
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics
include("required_scripts_2.jl")

## Specify the model
model = @reaction_network begin
   k₁,             ∅ --> A
    1.0-1.0*A/I,    A --> ∅
    k₂*A*A,         ∅ --> I
    k₂,             I --> ∅
end

## Specify model parameters
 params = model_parameters()
 N_screen = 120
params.reaction["k₁"] = screen_values(min = 0.01, max = 0.4, mode = "linear", number = N_screen) 
params.reaction["k₂"] = screen_values(min = 0.1, max = 2.0, mode = "linear", number = N_screen)   

params.diffusion["A"] = [0.1] 
params.diffusion["I"] = [1.0]

## Compile functions and solvers
include("simulate_scripts_2.jl")


## Find Turing params

@time turing_params = returnTuringParams(model, params);

## Plotting
using Plots
U₀ = VectorOfArray(turing_params.steady_state_values)[1,:]
V₀ = VectorOfArray(turing_params.steady_state_values)[2,:]
a = get_param(model, turing_params,"k₁","reaction")
b = get_param(model, turing_params,"k₂","reaction")
Dᵥ = get_param(model, turing_params,"I","diffusion")

#Define diffusion as a variable
d=params.diffusion["A"]
ã = range(0.01,stop=0.4,length=100)
Ũ₀ = ã / 5
Ṽ₀ = 1 .+ã.^2 ./ 25.
D̃ᵥ = Dᵥ[1]
i = findall(Dᵥ .== D̃ᵥ)
b̃H = (1 .-ã) ./ (1 .+ã)
b̃T = 2 ./d .+((1 .-ã)./(d .*(1 .+ã))) .-(2 ./d).*(2 ./(1 .+ã)).^0.5

Plots.plot(ã,[b̃H,b̃T], labels=["lower bound" "upper bound"], xlabel="σ", ylabel="k",linewidth = 5)
Plots.scatter!(a[i],b[i],strokewidth = 0, label = "numerical sol.", color = :black, markersize=1)
title!("Analytical and numerical solution for G-M model")
