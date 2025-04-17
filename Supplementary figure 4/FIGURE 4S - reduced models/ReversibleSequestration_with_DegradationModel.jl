# Figure 4S - Reversible Sequestration Model with degradation - showing that if deg(C)/k- is small enough, hIS tends to be negative be negative.
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
    k₋/Kₘ,                     A + I --> C
    k₋,                     C  --> A + I
    k₄,                     C --> ∅
    HillFunction(S,1,1,n₁),            ∅ --> A
    HillFunction(S,k₅,K₂,n₂),           ∅ --> I
    k₆*A,                           ∅ --> S 
end

## Specify model parameters
params = model_parameters()
N_screen = 4 
params.reaction["k₂"] = screen_values(min = .1, max = 10, mode = "log", number = 2)
params.reaction["k₃"] = [0.05]
params.reaction["k₄"] = screen_values(min = .01, max = 10, mode = "log", number = 4) 
params.reaction["k₅"] =screen_values(min = 1, max = 15, mode = "log", number = 2)
params.reaction["k₆"] = [0.1,10,100]
params.reaction["Kₘ"] =screen_values(min = 0.01, max = 100, mode = "log", number = 4) 
params.reaction["k₋"] =screen_values(min = .1, max = 100, mode = "log", number = 4) 
params.reaction["K₂"] = screen_values(min = .005, max = 5, mode = "log", number = 3)
params.reaction["n₁"] = screen_values(min = -3, max = 3, mode = "linear", number = 30)                                          
params.reaction["n₂"] =  screen_values(min = -3, max = 3, mode = "linear", number = 30)           

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = .01, max = 100, mode = "log", number = 3)
params.diffusion["C"] =  screen_values(min = .01, max = 100, mode = "log", number = 3)

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
Kₘ = get_param(model, turing_params,"Kₘ","reaction")
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
Plots.scatter(hAS, hIS, markersize=3, color=darkblue, alpha=1.0, markerstrokecolor=darkblue, legend=false, 
              xlims=(-3,3), ylims=(-3,3), 
              framestyle=:origin, tickfontsize=12, grid=false)
x1=[-5;1]
y1=[-1;-1]

x2=[1;1]
y2=[-1;-5]
lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 
Plots.plot!(x1,y1, linewidth = 4, color=lightblue)
Plots.plot!(x2,y2, linewidth = 4, color=lightblue)

x1=[1;5] 
y1 = (x1.-1)
x2=[1;1]
y2=[0;5]

Plots.plot!(x1,y1, linewidth = 4, color=lightblue)
Plots.plot!(x2,y2, linewidth = 4, color=lightblue)


##  test importance of kc/k-

ratio = k₄ ./ k₋
unique_values = unique(ratio)
println(unique_values)

subset1 = findall(ratio .== 0.00010000000000000002)
findall(hIS[subset1] .> 0.0)
findall(hIS[subset1] .< 0.0)


subset2 = findall((ratio .== 0.0010000000000000002) .| (ratio .== 0.001))
findall(hIS[subset2] .> 0.0)
findall(hIS[subset2] .< 0.0)

subset3 = findall((ratio .== 0.010000000000000002) .| (ratio .== 0.01))
findall(hIS[subset3] .> 0.0)
findall(hIS[subset3] .< 0.0)

subset4 = findall((ratio .== 0.10000000000000002) .| (ratio .== 0.1))
findall(hIS[subset4] .> 0.0)
findall(hIS[subset4] .< 0.0)

subset5 = findall(ratio .== 1.0)
findall(hIS[subset5] .> 0.0)
findall(hIS[subset5] .< 0.0)

subset6 = findall(ratio .== 10.0)
findall(hIS[subset6] .> 0.0)
findall(hIS[subset6] .< 0.0)

subset7 = findall(ratio .== 100.0)
findall(hIS[subset7] .> 0.0)
findall(hIS[subset7] .< 0.0)
