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
    HillFunction(S24n,1,1,n₁),            ∅ --> A 
    HillFunction(S24n,k₅,K₂,n₂),           ∅ --> I 
   A/I-kR*R,                                               ∅ --> R  
  730.0-kphos*R*S2c-kin2*S2c+kex2*S2n-0.1*S2c, ∅ --> S2c   
  kphos*R*S2c-kin2*pS2c-kon*pS2c*(S2c+2*pS2c)+koff*(S24c+2*S22c) +kex2*pS2n-0.1*pS2c, ∅ --> pS2c
  730.0-kin4*S4c-kon*pS2c*S4c+koff*S24c+kex4*S4n-0.1*S4c, ∅ --> S4c   
  kon*pS2c*S4c-koff*S24c-kin2*CIF*S24c-0.1*S24c, ∅ --> S24c
  kon*pS2c^2-koff*S22c-kin2*CIF*S22c-0.1*S22c, ∅ --> S22c   
  a*kin2*S2c-a*kex2*S2n+kdephos*PPase*pS2n -0.1*S2n, ∅ --> S2n    
  a*kin2*pS2c-a*kex2*pS2n-kdephos*PPase*pS2n-kon*pS2n*(S4n+2*pS2n)+koff*(S24n+2*S22n)-0.1*pS2n, ∅ --> pS2n
  a*kin4*S4c-a*kex4*S4n-kon*pS2n*S4n+koff*S24n-0.1*S4n, ∅ --> S4n 
  a*kin2*CIF*S24c+kon*pS2n*S4n-koff*S24n-0.1*S24n, ∅ --> S24n 
  a*kin2*CIF*S22c+kon*pS2n^2-koff*S22n-0.1*S22n, ∅ --> S22n    
end


## Specify model parameters
params = model_parameters()
params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.reaction["k₅"] =screen_values(min = 0.03, max = 30, mode = "log", number = 4)
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)                                 
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)      

params.reaction["kR"] = screen_values(min = 0.1, max = 10, mode = "log", number = 4)  

#other parameters
Tfactor = 1000
Mfactor = 0.01

params.reaction["kphos"] = [4.0*10^(-4)*Tfactor*Mfactor]
params.reaction["kdephos"] = [6.6*10^(-3)*Tfactor*Mfactor]
params.reaction["kin2"] = [2.6*10^(-3)*Tfactor]
params.reaction["kex2"] = [5.6*10^(-3)*Tfactor]
params.reaction["kin4"] = [2.6*10^(-3)*Tfactor]
params.reaction["kex4"] = [2.6*10^(-3)*Tfactor]
params.reaction["CIF"] = [5.7*1]   
params.reaction["kon"] = [1.8*10^(-3)*Tfactor*Mfactor]
params.reaction["koff"] = [1.6*10^(-2)*Tfactor]
params.reaction["a"] = [2.3*1]   
params.reaction["PPase"] = [1*1]

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.1, max = 10, mode = "log", number = 5) 


## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_general.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  

##
#Testing the necessary condition

@load "turing_paramsNS_wholePathway.jld2"
#First test:
A = StructArray((turing_params.reaction_params, turing_params.steady_state_values, turing_params.diffusion_constants))
B = StructArray((turing_params_NS.reaction_params, turing_params_NS.steady_state_values, turing_params_NS.diffusion_constants))

issubset(B,A)
setdiff(A,B)
