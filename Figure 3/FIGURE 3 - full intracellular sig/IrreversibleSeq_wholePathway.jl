# Code of Figure 3
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_2.jl")

include("plotting_scripts.jl")


function myOwnFunction(SIG,k,K,n)
    ratio = max((SIG),0)
    
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end

model = @reaction_network begin
    1,                      A --> ∅
    k₂,                     I --> ∅
    myOwnFunction(S14n,1,1,n₁) - k₊*A*I ,            ∅ --> A 
    myOwnFunction(S14n,k₅,K₂,n₂) - k₊*A*I,           ∅ --> I 
   A-kR*R,                                               ∅ --> R  
  730.0-kphos*R*S1c-kin2*S1c+kex2*S1n-0.1*S1c, ∅ --> S1c   #Constant is 73.0/Mfactor*degradation_rate
  kphos*R*S1c-kin2*pS1c-kon*pS1c*(S1c+2*pS1c)+koff*(S14c+2*S12c) +kex2*pS1n-0.1*pS1c, ∅ --> pS1c
  730.0-kin4*S4c-kon*pS1c*S4c+koff*S14c+kex4*S4n-0.1*S4c, ∅ --> S4c   
  kon*pS1c*S4c-koff*S14c-kin2*CIF*S14c-0.1*S14c, ∅ --> S14c
  kon*pS1c^2-koff*S12c-kin2*CIF*S12c-0.1*S12c, ∅ --> S12c   
  a*kin2*S1c-a*kex2*S1n+kdephos*PPase*pS1n -0.1*S1n, ∅ --> S1n    
  a*kin2*pS1c-a*kex2*pS1n-kdephos*PPase*pS1n-kon*pS1n*(S4n+2*pS1n)+koff*(S14n+2*S12n)-0.1*pS1n, ∅ --> pS1n
  a*kin4*S4c-a*kex4*S4n-kon*pS1n*S4n+koff*S14n-0.1*S4n, ∅ --> S4n 
  a*kin2*CIF*S14c+kon*pS1n*S4n-koff*S14n-0.1*S14n, ∅ --> S14n 
  a*kin2*CIF*S12c+kon*pS1n^2-koff*S12n-0.1*S12n, ∅ --> S12n    
end


## Specify model parameters
params = model_parameters()
N_screen = 4

params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3) 
params.reaction["k₅"] = screen_values(min = 1, max = 10, mode = "log", number = 2)  
params.reaction["k₊"] =screen_values(min = 1, max = 100, mode = "log", number = 6) 
params.reaction["K₂"] = screen_values(min = .1, max = 2, mode = "log", number = 3) 
params.reaction["n₁"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0]                                       
params.reaction["n₂"] = [0.5,1.0,1.3,1.6,2.0,2.3,2.6,3.0]     

params.reaction["kR"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3) 

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
params.diffusion["I"] =  screen_values(min = 0.3, max = 30, mode = "log", number = 3) 


## Compile functions and solvers
include("simulate_scripts_2.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  

## Compute hIS/hAS

k₂ = get_param(model, turing_params,"k₂","reaction")
k₅ = get_param(model, turing_params,"k₅","reaction")
k₊ = get_param(model, turing_params,"k₊","reaction")
K₂ = get_param(model, turing_params,"K₂","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")
A = VectorOfArray(turing_params.steady_state_values)[1,:]
I = VectorOfArray(turing_params.steady_state_values)[2,:]
S14n = VectorOfArray(turing_params.steady_state_values)[12,:]
kon = get_param(model, turing_params,"kon","reaction")
kphos = get_param(model, turing_params,"kphos","reaction")
S1c = VectorOfArray(turing_params.steady_state_values)[3,:]
pS1c = VectorOfArray(turing_params.steady_state_values)[4,:]

hAS = n₁ ./ (1 .+ S14n.^n₁)
hIS = n₂ .*K₂.^n₂ ./ (K₂.^n₂ .+ S14n.^n₂)

darkblue= RGBA(0/255, 114/255, 178/255, 1.0)
Plots.scatter(hAS, hIS, markersize=4, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-0.3,3), ylims=(-0.3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)

x1=[1;5] 
y1 = (x1.-1)
x2=[1;1]
y2=[0;5]

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1, y1, linewidth=9, color=lightblue)
Plots.plot!(x2,y2, linewidth = 9, color=lightblue)





