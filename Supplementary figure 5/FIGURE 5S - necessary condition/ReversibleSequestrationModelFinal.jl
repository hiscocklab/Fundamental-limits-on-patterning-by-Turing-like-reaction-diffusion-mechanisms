## Figure 5S - necessary condition
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_NECESSARY_CONDITION.jl")

include("plotting_scripts.jl")

function HillFunction(S,k,K,n)
    ratio = max((S),0)
    
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end

model = @reaction_network begin
    1,                     A --> ∅
    k₂,                     I --> ∅
    k₊*A*I - k₋*C,                              ∅ --> C
    k₃,                     S --> ∅
    HillFunction(S,1,1,n₁) - k₊*A*I + k₋*C,            ∅ --> A
    HillFunction(S,k₅,K₂,n₂)- k₊*A*I + k₋*C,           ∅ --> I
    k₆*A,                           ∅ --> S 
end

## Specify model parameters
params = model_parameters()
N_screen = 4
params.reaction["k₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)
params.reaction["k₃"] = [0.01]
params.reaction["k₅"] =screen_values(min = 0.1, max = 10, mode = "log", number = 5)
params.reaction["k₆"] = screen_values(min = .01, max = 10, mode = "log", number = 5)
params.reaction["k₊"] =screen_values(min = 0.1, max = 100, mode = "log", number = 5)
params.reaction["k₋"] = screen_values(min = .1, max = 100, mode = "log", number = 5)
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5)
params.reaction["n₁"] = screen_values(min = -3, max = -0.1, mode = "linear", number = 10)                                          
params.reaction["n₂"] =  screen_values(min = -3, max = -1, mode = "linear", number =8)           

params.diffusion["A"] = [1.0] #single values must be in square brackets
params.diffusion["I"] =  screen_values(min = .1, max = 10, mode = "log", number = 5)
params.diffusion["C"] =  screen_values(min = .1, max = 10, mode = "log", number = 5)
## Compile functions and solvers
include("simulate_scripts_NECESSARY_CONDITION_4_BY_4.jl") 

## Turing screen
@time turing_params = returnTuringParams(model, params);

## Compute Robustness for all varying parameters
model.ps
model.states
params.diffusion
using Measures
ps,ds,ics, _, _, _ = returnParameterSets(model, params)
n_iter = prod([length.(ps); length.(ics); length.(ds)])
i = ["k₂", "k₊", "k₋", "k₃", "n₁", "k₅", "K₂", "n₂", "k₆"]  #Order must match order of turing_params.reaction_params (model.ps)

Robustness= zeros(length(i), 0) 
 # Initialize with 0 columns
for idx_i in 1:length(i)
    j = params.reaction[i[idx_i]]
    
    # Update the number of columns in Dadada
    Robustness = zeros(1, size(Robustness, 2) + length(j))
    
    for idx_j in 1:length(j)
        # Construct the dynamic indices
        reaction_name = i[idx_i]

        # Count occurrences and update Robustness
        Robustness[1, end - length(j) + idx_j] = round(100*count(x -> x == params.reaction[reaction_name][idx_j], VectorOfArray(turing_params.reaction_params)[idx_i, :]) * length(params.reaction[reaction_name]) / n_iter;digits = 1)
    end
    
    p = Plots.plot(
        Plots.heatmap(
            reshape(Robustness, 1, length(Robustness)),
            color=:Reds,
            cmin=0,
            cmax=1,
            xlabel=i[idx_i],
            title="Robustness Heatmap",
            yflip=true,
            yticks=false,
            xticks=(1:length(j), [string(round(val, digits=2)) for val in j]),  # Format tick labels to two decimal places
            size=(500, 160),  # Adjust the overall size of the plot
            legend=:top,
            tickfontsize=8,
            guidefontsize=10,
            bottom_margin=10Measures.mm,
            right_margin=10Measures.mm,
            clim = (0,1.0)
        )
    )
display(p)
Robustness= zeros(length(i), 0)
end 
#
#Robustness for diffusion: NOTICE: Order i the same way as model.states is ordered
turing_params.diffusion_constants
#turing_params.diffusion_constants
#params.diffusion # Probelsm? Order!

i = ["A";"I";"C";"S"]
Robustness= zeros(length(i), 0) 
 # Initialize with 0 columns
for idx_i in 1:length(i)
    j = params.diffusion[i[idx_i]]
    
    # Update the number of columns in Dadada
    Robustness = zeros(1, size(Robustness, 2) + length(j))
    
    for idx_j in 1:length(j)
        # Construct the dynamic indices
        reaction_name = i[idx_i]

        # Count occurrences and update Robustness
        Robustness[1, end - length(j) + idx_j] = round(100*count(x -> x == params.diffusion[reaction_name][idx_j], VectorOfArray(turing_params.diffusion_constants)[idx_i, :]) * length(params.diffusion[reaction_name]) / n_iter;digits = 1)
    end
    
    p = Plots.plot(
        Plots.heatmap(
            reshape(Robustness, 1, length(Robustness)),
            color=:Reds,
            cmin=0,
            cmax=1,
            xlabel=i[idx_i],
            title="Robustness Heatmap",
            yflip=true,
            yticks=false,
            xticks=(1:length(j), [string(round(val, digits=2)) for val in j]),  # Format tick labels to two decimal places
            size=(500, 160),  # Adjust the overall size of the plot
            legend=:top,
            tickfontsize=8,
            guidefontsize=10,
            bottom_margin=10Measures.mm,
            right_margin=10Measures.mm,
            clim=(0,1)
        )
    )
display(p)
Robustness= zeros(length(i), 0)
end 
