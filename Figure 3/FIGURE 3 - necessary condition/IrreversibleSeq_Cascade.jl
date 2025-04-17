# Code of Figure 3 - necessary condition
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
    0.01,                      S9 --> ∅    
    1,                     a --> ∅
    1,                     i --> ∅
    1,                     S1 --> ∅
    1,                     S2 --> ∅
    1,                     S3 --> ∅
    1,                     S4 --> ∅
    1,                     S5 --> ∅
    1,                     S6 --> ∅
    1,                     S7 --> ∅
    1,                     S8 --> ∅
    kR,                    R --> ∅
    a - k₊*A*I ,            ∅ --> A  
    i - k₊*A*I,           ∅ --> I 
    A,                           ∅ --> R 
    1*R,                           ∅ --> S1
    1*S1,                           ∅ --> S2
    1*S2,                           ∅ --> S3
    1*S3,                           ∅ --> S4
    1*S4,                           ∅ --> S5
    1*S5,                           ∅ --> S6  
    1*S6,                           ∅ --> S7
    1*S7,                           ∅ --> S8
    k₆*S8,                            ∅ --> S9
    HillFunction(S9,1,1,n₁),         ∅ --> a
    HillFunction(S9,k₅,K₂,n₂),        ∅ --> i
end


## Specify model parameters
params = model_parameters()

params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.reaction["k₅"] = screen_values(min = 1, max = 10, mode = "log", number = 2)  
params.reaction["k₆"] = [3.0]
params.reaction["k₊"] = screen_values(min = 1, max = 100, mode = "log", number = 6)
params.reaction["K₂"] = screen_values(min = .1, max = 2, mode = "log", number = 3)
params.reaction["n₁"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0]                             
params.reaction["n₂"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0] 
params.reaction["kR"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.3, max = 30, mode = "log", number = 3)  

# k₃ chosen to be 0.01

## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  
