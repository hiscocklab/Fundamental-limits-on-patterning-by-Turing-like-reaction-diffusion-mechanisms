# Code for Figure 3S - necessary condition
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
    0.1,                      S9 --> ∅    
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
    a ,            ∅ --> A  
    i ,           ∅ --> I 
    A/I,                           ∅ --> R 
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
params.reaction["k₅"] =screen_values(min = 0.03, max = 30, mode = "log", number = 4)
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)                                 
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)  
params.reaction["k₆"] =[0.7]
params.reaction["kR"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.1, max = 10, mode = "log", number = 5)  

# k₃ chosen to be 0.01

## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")   

## Turing screen
@time turing_params = returnTuringParams(model, params);  

##
#Testing if every necessary condition satisfying set is in Turing-pattern-forming set

@load "turing_paramsNS_Cascade.jld2"
#First test:
A = StructArray((turing_params.reaction_params, turing_params.steady_state_values, turing_params.diffusion_constants))
B = StructArray((turing_params_NS.reaction_params, turing_params_NS.steady_state_values, turing_params_NS.diffusion_constants))

issubset(B,A)
setdiff(A,B)