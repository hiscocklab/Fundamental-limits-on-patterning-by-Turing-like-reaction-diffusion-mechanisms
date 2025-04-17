# Code for Figure 3 - necessary condition
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
    1,                      A --> ∅
    k₂,                     I --> ∅
    k₃,                      S --> ∅
    kR,                     R --> ∅
    1,                     a --> ∅
    1,                     i --> ∅
    a - k₊*A*I ,            ∅ --> A  
    i - k₊*A*I,           ∅ --> I 
    A,                           ∅ --> R 
    k₆*R,                            ∅ --> S
    HillFunction(S,1,1,n₁),         ∅ --> a
    HillFunction(S,k₅,K₂,n₂),        ∅ --> i
end


## Specify model parameters
params = model_parameters()

params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3) 
params.reaction["k₅"] = screen_values(min = 1, max = 10, mode = "log", number = 2)   
params.reaction["k₆"] = [0.4]
params.reaction["k₊"] = screen_values(min = 1, max = 100, mode = "log", number = 6) 
params.reaction["K₂"] = screen_values(min = .1, max = 2, mode = "log", number = 3) 
params.reaction["n₁"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0]                                 
params.reaction["n₂"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0]
params.reaction["kR"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.3, max = 30, mode = "log", number = 3) 

#Varying across 
params.reaction["k₃"] = screen_values(min = 0.01, max = 10, mode = "log", number = 4) #when running I had to split this in four runs

## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  
