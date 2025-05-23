## Figure 4S - necessary condition
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_NECESSARY_CONDITION.jl")
include("plotting_scripts.jl")

function HillFunction(SIG,k,K,n)
    ratio = max((SIG),0)
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end
model = @reaction_network begin
    1,                     A1 --> ∅
    k₂,                     A2 --> ∅
    k₃,                     SA1 --> ∅
    k₄,                     SA2 --> ∅
    HillFunction(SA1,1,1,n₁)*myOwnFunction(SA2,1,1,n₂),            ∅ --> A1
    HillFunction(SA1,k₅,K₃,n₃),           ∅ --> A2
    k₆*A1,                           ∅ --> SA1
    k₇*A2,                              ∅ --> SA2
end

## Specify model parameters
params = model_parameters()
N_screen = 5
params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 4)
params.reaction["k₃"] = screen_values(min = .01, max = 10, mode = "log", number = 4)
params.reaction["k₄"] = screen_values(min = 0.01, max = 10, mode = "log", number = 4)
params.reaction["k₅"] =screen_values(min = 0.1, max = 10, mode = "log", number = 4)
params.reaction["k₆"] = screen_values(min = 0.3, max = 7, mode = "log", number = 3)
params.reaction["k₇"] = screen_values(min = 0.3, max = 10, mode = "log", number = 3)
params.reaction["K₃"] = screen_values(min = 0.05, max = 2, mode = "log", number = 4)
params.reaction["n₁"] = [1.5]
params.reaction["n₂"] =  screen_values(min = -3, max = 3, mode = "linear", number = 13)
params.reaction["n₃"] =  screen_values(min = -3, max = 3, mode = "linear", number = 13 )

params.diffusion["A1"] = [1.0] 
params.diffusion["A2"] =  screen_values(min = .01, max = 70, mode = "log", number = 5)

## Compile functions and solvers 
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen  
@time turing_params = returnTuringParams(model, params);
