## import and load
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics
include("required_scripts_2.jl")

## Specify the model
model = @reaction_network begin
    γ*a + γ*U^2*V,  ∅ --> U
    γ,              U --> ∅
    γ*b,            ∅ --> V
    γ*U^2,          V --> ∅
end

## Specify model parameters
params = model_parameters()
N_screen = 100
params.reaction["a"] = screen_values(min = 0, max = 0.5, mode = "linear", number = N_screen)
params.reaction["b"] = screen_values(min = 0, max = 3, mode = "linear", number = N_screen)
params.reaction["γ"] = [1]
params.diffusion["U"] = [1.0] 
params.diffusion["V"] = [50.0]

## Compile functions and solvers
include("simulate_scripts_2.jl")

## Find Turing params
@time turing_params = returnTuringParams(model, params);

##Plotting 

using Plots

a = get_param(model, turing_params,"a","reaction")
b = get_param(model, turing_params,"b","reaction")
b̃ = range(0,stop=3,length=100)
d = params.diffusion["V"]/params.diffusion["U"]

b̃H = (1 ./ (1.44 .* (-9*b̃ .+ (sqrt.(3 .+ 81*b̃.^2))).^(1/3))) .- b̃ .- ((-9*b̃ .+ (sqrt.(3 .+ 81*b̃.^2))).^(1/3)) ./ (2.08)
h = (27*b̃.*d .+ d.^(3/2) .+ sqrt.(729*b̃.^2 .* d.^2 .+ 54*b̃.*d.^(5/2))).^(1/3)
b̃T = (1/3).*(-3*b̃ .- 2*sqrt.(d) .+ d./h .+ h )

Plots.plot(b̃,[b̃H,b̃T], labels=["lower bound" "upper bound"], ylimits =(0,0.45), xlimits =(0, 3), xlabel="b", ylabel="a", linewidth = 5)
Plots.scatter!(b,a,strokewidth = 0,label = "numerical sol.", color = :black, markersize=1)
title!("Analytical and numerical sol. - Schnakenberg model")

