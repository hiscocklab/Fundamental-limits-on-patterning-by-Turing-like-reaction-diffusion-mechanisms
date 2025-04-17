## Figure 4

using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_2.jl")
include("plotting_scripts.jl")

function HillFunction(S,k,K,n)
    ratio = abs(S)
    
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end

model = @reaction_network begin
    1,                     A --> ∅
    k₂,                     I --> ∅
    k₃,                     S --> ∅
    HillFunction(S,1,1,n₁),        ∅ --> A
    HillFunction(S,k₅,K₂,n₂),      ∅ --> I
    (k₆*(A/I)),                 ∅ --> S 
end

## Specify model parameters
params = model_parameters()
N_screen = 5
params.reaction["k₂"] = screen_values(min = .1, max = 10, mode = "log", number = N_screen)
params.reaction["k₃"] = screen_values(min = .01, max = 10, mode = "log", number = 4)
params.reaction["k₅"] = screen_values(min = .1, max = 10, mode = "log", number = N_screen)
params.reaction["k₆"] = screen_values(min = .001, max = 3, mode = "log", number = N_screen)
params.reaction["K₂"] = screen_values(min = .1, max = 5, mode = "log", number = N_screen)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 21)                                        
params.reaction["n₂"] = screen_values(min = -3, max = 3, mode = "linear", number = 21)           

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = .01, max = 100, mode = "log", number = 5)

## Compile functions and solvers 
include("simulate_scripts_2.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);

## Compute and Plot hIS/hAS

k₂ = get_param(model, turing_params,"k₂","reaction")
k₃ = get_param(model, turing_params,"k₃","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
k₅ = get_param(model, turing_params,"k₅","reaction")
K₂ = get_param(model, turing_params,"K₂","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")
k₆ = get_param(model, turing_params,"k₆","reaction")

#diffusion
D_I =get_param(model, turing_params,"I","diffusion")
# steady states
A = VectorOfArray(turing_params.steady_state_values)[1,:]
I = VectorOfArray(turing_params.steady_state_values)[2,:]
S = VectorOfArray(turing_params.steady_state_values)[3,:]

hAS = n₁ ./ (1 .+ S .^n₁)
hIS =  K₂ .^ n₂ .* n₂ ./(K₂.^n₂ .+ S.^n₂) 

darkblue= RGBA(0/255, 114/255, 178/255, 1.0)
Plots.scatter(hAS, hIS, markersize=2, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-3,3), ylims=(-3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)
#necessary conditions
x=[-7;-3;-2;-1;0;1;2;3;7]
x1_index = x .<1
x1 =x[x1_index]
x2_index = x .>=1
x2 =x[x2_index]

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1,(x1 .- 1), linewidth = 6, color=lightblue)
Plots.plot!(x2,x2 .- 1, linewidth = 6, color=lightblue)
Plots.plot!(x1,-ones(length(x1)),linewidth=6, color=lightblue) 
Plots.plot!(ones(length(x2)),x2 .-1, linewidth=6, color=lightblue)
