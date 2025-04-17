## import and load
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics
include("required_scripts_2.jl")

## Specify the model
model = @reaction_network begin
    (k₂*1.0+0.1*C*C*D),  ∅ --> A
    A*A*B,             ∅ --> A
    1.0,               A --> ∅
    k₁,                ∅ --> B
    A*A,               B --> ∅
    (k₂*1.0+0.1*A*A*B),  ∅ --> C
    C*C*D,             ∅ --> C
    1.0,               C --> ∅
    k₁,                ∅ --> D
    C*C,               D --> ∅
end

## Specify model parameters
params = model_parameters()
Nscreen = 100

params.reaction["k₁"] = screen_values(min = 0.001, max = 4.0, mode = "linear", number = N_screen)
params.reaction["k₂"] = screen_values(min = 0.001, max = 1.0, mode = "linear", number = N_screen)

params.diffusion["A"] = [1.0] 
params.diffusion["B"] = [50.0]
params.diffusion["C"] = [1.0] 
params.diffusion["D"] = [50.0]


## Compile functions and solvers
include("simulate_scripts_2.jl")

## Test for Turing

@time turing_params = returnTuringParams(model, params);

## Plotting

using Plots
U₀ = VectorOfArray(turing_params.steady_state_values)[1,:]
V₀ = VectorOfArray(turing_params.steady_state_values)[2,:]
a = get_param(model, turing_params,"k₁","reaction")
b = get_param(model, turing_params,"k₂","reaction")
Dᵥ = get_param(model, turing_params,"B","diffusion")

b̃ = range(0.001, stop = 4, length = 100)
ã = range(0.01,stop=0.5,length=30)
Ũ₀ = ã / 5
Ṽ₀ = 1 .+ã.^2 ./ 25.
D̃ᵥ = Dᵥ[1]
d=50
i = findall(Dᵥ .== D̃ᵥ)
h = (27*b̃ .+ 3*sqrt.(81*b̃.^2 .+3)).^(1/3)
b̃H2 = (3 .^(1/3).*10 .^(2/3).*(99*b̃ .+ 3 .^(1/2).*(3267*b̃.^2 .+ 100).^(1/2)).^(1/3))./30 .- (11*b̃)./10 .- (3 .^(2/3).*10 .^(1/3))./(3 .*(99*b̃ .+ 3 .^(1/2).*(3267*b̃.^2 .+ 100).^(1/2)).^(1/3))
b̃T3 = (10 .^(2/3).*(297*b̃.*d .+ 10*d.^(3/2) .+ 3*33 .^(1/2).*(b̃*d.^2 .*(297*b̃ .+ 20*d.^(1/2))).^(1/2)).^(1/3))./30 .- (11*b̃)./10 .- (2*d.^(1/2))./3 .+ (10 .^(1/3).*d)./(3 .*(297*b̃.*d .+ 10*d.^(3/2) .+ 3*33 .^(1/2).*(b̃*d.^2 .*(297*b̃ .+ 20*d.^(1/2))).^(1/2)).^(1/3))

Plots.plot(b̃,[b̃H2,b̃T3],labels=["lower bound" "upper bound"],ylimits =(0,1), xlimits =(0, 3),xlabel="b", ylabel="a",linewidth = 5)

Plots.scatter!(a[i],b[i],strokewidth = 0, label = "numerical sol.", color = :black, markersize=1)
title!("Analytical and numerical sol. - C1 model")

