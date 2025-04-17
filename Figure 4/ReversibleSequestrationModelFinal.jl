## Figure 4

using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_2.jl")

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
    HillFunction(S,1,1,n₁) - k₊*A*I + k₋*C,            ∅ --> A
    HillFunction(S,k₅,K₂,n₂)- k₊*A*I + k₋*C,           ∅ --> I
    k₆*A,                           ∅ --> S 
    k₊*A*I - k₋*C,                              ∅ --> C
end

## Specify model parameters
params = model_parameters()

params.reaction["k₂"] = screen_values(min = .1, max = 10, mode = "log", number = 2)
params.reaction["k₃"] = screen_values(min = .01, max = 10, mode = "log", number = 3)
params.reaction["k₅"] =screen_values(min = 5, max = 15, mode = "log", number = 2)
params.reaction["k₆"] = [0.1,10,100] 
params.reaction["k₊"] =screen_values(min = 10, max = 100, mode = "log", number = 2)
params.reaction["k₋"] = screen_values(min = .1, max = 10, mode = "log", number = 2)
params.reaction["K₂"] = screen_values(min = .001, max = 2, mode = "log", number = 4)  
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 50)                                          
params.reaction["n₂"] =  screen_values(min = -3, max = 3, mode = "linear", number = 50)           
 
params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = .01, max = 100, mode = "log", number = 4)
params.diffusion["C"] =  screen_values(min = .01, max = 100, mode = "log", number = 4)

## Compile functions and solvers  
include("simulate_scripts_2.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);

## Compute and Plot hIS/hAS
k₂ = get_param(model, turing_params,"k₂","reaction")
k₃ = get_param(model, turing_params,"k₃","reaction")
k₅ = get_param(model, turing_params,"k₅","reaction")
k₆ = get_param(model, turing_params,"k₆","reaction")
k₊ = get_param(model, turing_params,"k₊","reaction")
k₋ = get_param(model, turing_params,"k₋","reaction")
K₂ = get_param(model, turing_params,"K₂","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")

A = VectorOfArray(turing_params.steady_state_values)[1,:]
I = VectorOfArray(turing_params.steady_state_values)[2,:]
S = VectorOfArray(turing_params.steady_state_values)[3,:]
C = VectorOfArray(turing_params.steady_state_values)[4,:]

hAS = n₁ ./ (1 .+ S.^n₁) 
hIS = n₂ .*K₂.^n₂ ./ (K₂.^n₂ .+ S.^n₂)


darkblue= RGBA(0/255, 114/255, 178/255, 1.0)
Plots.scatter(hAS, hIS, markersize=2, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-3,3), ylims=(-3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)
# Plot necessary condition:
x1=[-5;1] 
y1 = [-1;-1]
x2=[1;1]
y2=[-1;-5]

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1, y1, linewidth=6, color=lightblue)
Plots.plot!(x2,y2, linewidth = 6, color=lightblue)