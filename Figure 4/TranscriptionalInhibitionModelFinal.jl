## Figure 4
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
params.reaction["n₃"] =  screen_values(min = -3, max = 3, mode = "linear", number = 13)

params.diffusion["A1"] = [1.0] 
params.diffusion["A2"] =  screen_values(min = .01, max = 70, mode = "log", number = 5)

## Compile functions and solvers
include("simulate_scripts_2.jl")

## Turing screen  
@time turing_params = returnTuringParams(model, params);

## Compute and Plot hIS/hAS
k₂ = get_param(model, turing_params,"k₂","reaction")
k₃ = get_param(model, turing_params,"k₃","reaction")
k₄ = get_param(model, turing_params,"k₄","reaction")
k₅ = get_param(model, turing_params,"k₅","reaction")
k₆ = get_param(model, turing_params,"k₆","reaction")
k₇ = get_param(model, turing_params,"k₇","reaction")
K₃ = get_param(model, turing_params,"K₃","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")
n₃ = get_param(model, turing_params,"n₃","reaction")

A1 = VectorOfArray(turing_params.steady_state_values)[1,:]
A2 = VectorOfArray(turing_params.steady_state_values)[2,:]
SA1 = VectorOfArray(turing_params.steady_state_values)[3,:]
SA2 = VectorOfArray(turing_params.steady_state_values)[4,:]

hASA2 = n₂ ./ (1 .+ SA2.^n₂)

hISA1 = n₃ .*K₃.^n₃ ./ (K₃.^n₃ .+ SA1.^n₃)

darkblue= RGBA(0/255, 114/255, 178/255, 1.0)
Plots.scatter(hASA2, hISA1, markersize=2, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-3,3), ylims=(-3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)
x1=[-5;5]
y1=[0;0]
x2=[0;0]
y2=[-5;5]

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1,y1, linewidth = 6, color=lightblue)
Plots.plot!(x2,y2, linewidth = 6, color=lightblue)
