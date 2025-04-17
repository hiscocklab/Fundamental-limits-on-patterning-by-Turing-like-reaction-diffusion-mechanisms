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
    1,                     A --> ∅
    k₂,                     I --> ∅
    k₃,                     S --> ∅
    HillFunction(S,1,1,n₁),            ∅ --> A
    HillFunction(S,k₅,K₂,n₂) ,           ∅ --> I
    k₆*A/I,                           ∅ --> S 
end


## Specify model parameters
params = model_parameters()
params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)  
params.reaction["k₅"] =screen_values(min = 0.03, max = 30, mode = "log", number = 4)  
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)   
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)                                 
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)      

params.reaction["k₆"] = screen_values(min = 0.01, max = 10, mode = "log", number = 6)  

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.1, max = 10, mode = "log", number = 5) 

#degradation of R and S
params.reaction["k₃"] = screen_values(min = 0.01, max = 10, mode = "log", number = 4) 

## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);



##
#Testing if necessary is sufficient (testing very precisely)

@load "turing_paramsNS.jld2" #saved Turing-pattern-forming parameter sets
A = StructArray((turing_params.reaction_params, turing_params.steady_state_values, turing_params.diffusion_constants))
B = StructArray((turing_params_NS.reaction_params, turing_params_NS.steady_state_values, turing_params_NS.diffusion_constants))

issubset(B,A)
setdiff(A,B)
