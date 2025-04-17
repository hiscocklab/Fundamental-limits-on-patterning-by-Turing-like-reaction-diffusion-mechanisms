## Figure 4S - necessary condition

using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_NECESSARY_CONDITION.jl")

include("plotting_scripts.jl")


function HillFunction(S,k,K,n)
    ratio = max((S),0)
    
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end
model = @reaction_network begin
    1,                     A --> ∅
    k₂,                     I --> ∅
    k₃,                     S --> ∅
    HillFunction(S,1,1,n₁) - k₊*A*I ,            ∅ --> A
    HillFunction(S,k₅,K₂,n₂)- k₊*A*I ,           ∅ --> I
    k₆*A,                           ∅ --> S 
end


## Specify model parameters
params = model_parameters()
N_screen = 4
params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.reaction["k₃"] = screen_values(min = .01, max = 10, mode = "log", number = 5)
params.reaction["k₅"] =screen_values(min = 1, max = 10, mode = "log", number = 3)
params.reaction["k₆"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)  
params.reaction["k₊"] =screen_values(min = 0.5, max = 100, mode = "log", number = 7)
params.reaction["K₂"] = screen_values(min = .1, max = 3, mode = "log", number = 3)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 21)                                        
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 21)          

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = .01, max = 100, mode = "log", number = 5)

## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);
