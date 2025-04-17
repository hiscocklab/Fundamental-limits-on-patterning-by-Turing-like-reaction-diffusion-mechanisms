# Code for Figure 3S
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_2.jl")

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
    a,            ∅ --> A  
    i,           ∅ --> I 
    A/I,                           ∅ --> R 
    k₆*R,                            ∅ --> S
    HillFunction(S,1,1,n₁),         ∅ --> a
    HillFunction(S,k₅,K₂,n₂),        ∅ --> i
end

## Specify model parameters
params = model_parameters()

params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 3)
params.reaction["k₅"] =screen_values(min = 0.03, max = 30, mode = "log", number = 4)
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)                                 
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 13)   
params.reaction["k₆"] =[0.07]
params.reaction["kR"] =screen_values(min = 0.1, max = 10.0 , mode = "log", number = 3)

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = 0.1, max = 10, mode = "log", number = 5) 

#Varying across 
params.reaction["k₃"] = screen_values(min = 0.1, max = 1.0, mode = "log", number = 2) 

## Compile functions and solvers
include("simulate_scripts_2.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  

## Compute and Plot hIS/hAS


k₂ = get_param(model, turing_params,"k₂","reaction")
k₃ = get_param(model, turing_params,"k₃","reaction")
k₅ = get_param(model, turing_params,"k₅","reaction")
k₆ = get_param(model, turing_params,"k₆","reaction")
K₂ = get_param(model, turing_params,"K₂","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")

A = VectorOfArray(turing_params.steady_state_values)[1,:]
I = VectorOfArray(turing_params.steady_state_values)[2,:]
S = VectorOfArray(turing_params.steady_state_values)[3,:]
R = VectorOfArray(turing_params.steady_state_values)[4,:]
a = VectorOfArray(turing_params.steady_state_values)[5,:]
i = VectorOfArray(turing_params.steady_state_values)[6,:]

hAS = n₁ ./ (1 .+ S.^n₁) 
hIS = n₂ .*K₂.^n₂ ./ (K₂.^n₂ .+ S.^n₂)


darkblue= RGBA(0/255, 114/255, 178/255, 1.0)
Plots.scatter(hAS, hIS, markersize=3, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-3,3), ylims=(-3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)
#Plot necessary condition:
x1=[1;5] 
y1 = (x1.-1)
x2=[1;1]
y2=[0;5]

x3=[-3;0] 
y3 = (x3.-1)
x4=[-3;0]
y4=[-1;-1]

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1, y1, linewidth=5, color=lightblue)
Plots.plot!(x2,y2, linewidth = 5, color=lightblue)
Plots.plot!(x3, y3, linewidth=5, color=lightblue)
Plots.plot!(x4,y4, linewidth = 5, color=lightblue)
