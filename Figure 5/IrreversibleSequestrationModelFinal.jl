# Figure 5
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load
include("required_scripts_2.jl")
include("plotting_scripts.jl")

function HillFunction(S,k,K,n)
    ratio = max(S,0)
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end
model = @reaction_network begin
    1,                     A --> ∅
    k₂,                     I --> ∅
    k₃,                     S --> ∅
    HillFunction(S,1,1,n₁) - k₊*A*I ,            ∅ --> A
    HillFunction(S,k₅,K₂,n₂)- k₊*A*I ,           ∅ --> I
    k₆*A,                           ∅ --> S 
end


## Specify model parameters
params = model_parameters()
N_screen = 5
params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 5) # deg I/deg A
params.reaction["k₃"] = [0.01]
params.reaction["k₅"] =screen_values(min = 0.1, max = 10, mode = "log", number = 5) # max transcription rate I/A
params.reaction["k₆"] = screen_values(min = 0.01, max = 10, mode = "log", number = 5) 
params.reaction["k₊"] =screen_values(min = 0.1, max = 10, mode = "log", number = 4)
params.reaction["K₂"] = screen_values(min = .1, max = 10, mode = "log", number = 5) #ligand concentration producing half occupation I/A
params.reaction["n₁"] = screen_values(min = 0, max = 3, mode = "linear", number = 10)                                        
params.reaction["n₂"] = screen_values(min = 0, max = 3, mode = "linear", number = 10)          

params.diffusion["A"] = [1.0] 
params.diffusion["I"] =  screen_values(min = .1, max = 10, mode = "log", number = 5)  # Diffusion constant I/A

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
K₂ = get_param(model, turing_params,"K₂","reaction")
n₁ = get_param(model, turing_params,"n₁","reaction")
n₂ = get_param(model, turing_params,"n₂","reaction")
A = VectorOfArray(turing_params.steady_state_values)[1,:]
I = VectorOfArray(turing_params.steady_state_values)[2,:]
S = VectorOfArray(turing_params.steady_state_values)[3,:]

hAS = n₁ ./ (1 .+ S.^n₁) # notice that hSA=1
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

lightblue = RGBA(86/255, 180/255, 233/255, 1.0) 

Plots.plot!(x1, y1, linewidth=6, color=lightblue)
Plots.plot!(x2,y2, linewidth = 6, color=lightblue)


## Compute Range Robustness
function returnAllVariables(model, params)
    ps = Array{Array{Float64,1},1}(undef,0)
    for p in parameters(model) 
        push!(ps,get!(params.reaction,string(p),[0])) ## default is zero
    end

    ds = Array{Array{Float64,1},1}(undef,0)
    ics = Array{Array{Float64,1},1}(undef,0)
    for state in states(model)  
        push!(ds,get!(params.diffusion,chop(string(state), head=0,tail=3),[0])) 
        push!(ics,get!(params.initial_condition,chop(string(state), head=0,tail=3),[1])) 
    end
    seeds = Int64.(params.random_seed)
    noise = params.initial_noise
    ls = params.domain_size
    return ps,ds,ics, ls, seeds, noise
end

function trimVariables(states)
    return chop.(string.(states), head=0,tail=3)
end
model.ps
model.states
function computeRobustness(params, turing_params)
    ps,ds,_, _, _, _ = returnAllVariables(model, params)
    full_ps = [ps; ds]
    full_labels = [string.(model.ps); string.("D_",trimVariables(model.states))]
    P = [Array(VectorOfArray(turing_params.reaction_params)); Array(VectorOfArray(turing_params.diffusion_constants))]
    
    #In general: get rid of n_ and all non-varying params
    P = P[[1,4,5,6,8,10],:]   
    full_labels = full_labels[[1,4,5,6,8,10]]
    full_ps = full_ps[[1,4,5,6,8,10]]

    robustness = zeros(Float64,size(P))
    
    N_p = size(P,1)

    Threads.@threads for i = 1:N_p
        Pⱼ = P[1:N_p .!= i,:]
        for i_n in 1:size(P)[2]
            pⱼ = Pⱼ[:,i_n]

            i_s = findall(i->all(j->pⱼ[j] == Pⱼ[j,i],1:size(Pⱼ,1)),1:size(Pⱼ,2))

            pmax = maximum(P[i,i_s]) 
            pmin = minimum(P[i,i_s])  
            
            robustness[i,i_n] = pmax/pmin
           
        end
    end
    return log10.(robustness), full_labels 
  
end
robustness, full_labels = computeRobustness(params, turing_params)
geometric_average_robustness = vec(10 .^(mean((log10.(robustness)),dims=1))) 

##
#Finding concrete parameter set to analyse

non_oscillating = findall(turing_params.non_oscillatory)   #just a check
subset1 = sortperm(geometric_average_robustness, rev=true)[1:10] # choose n most robust param sets
maximum(geometric_average_robustness)

subset2 = findall(geometric_average_robustness .>0) # above a certain robustness
subset3 = findall(string.(turing_params.pattern_phase) .== string([1,1,1])) # of a particular phase, look at unique(turing_params.phase)

subset = intersect(subset1)

index = subset[9]   #This turing parameter set will be analysed
param1 = get_params(model, turing_params[index])  #to agree with the robustness barcode
param1.domain_size = 45
sol = simulate(model,param1)

## Plot solution
plot_pattern(model,sol)

## Plot mRNA (add mRNA in the plot)

S = last(sol)[:,3]     #Check from model.states
x = range(0,1, length = size(S,1))
mRNA_A = 1 .* S .^n₁[index] ./(S .^n₁[index] .+ 1)
mRNA_A = mRNA_A ./ maximum(mRNA_A,dims=1)
mRNA_I =k₅[index] .* S .^ n₂[index] ./ (S .^ n₁[index] .+ K₂[index] .^n₂[index])
mRNA_I = mRNA_I ./ maximum(mRNA_I,dims=1)
plot!(x,mRNA_A,label = "mRNA_A", ticks=:none, linewidth=2)
plot!(x,mRNA_I,label = "mRNA_I", ticks=:none,linewidth=3)


## Plot dynamics  
plot_dynamics(model,sol)
 

##Save solution for a later date
@save "LigandSequestrationModelDEG_ONLY1RD.jld2" model params turing_params sol

## load saved solution and explore

@load "LigandSequestrationModelDEG_ONLY1RD.jld2"



##Robustness barcode for the concrete solution
model.states

function plot_robustness_barcode(turing_params, params, model, index) 
    custom_labels = ["k₂";"k₃";"n₁";"k₊"; "k₅"; "K₂"; "n₂";"k₆"] 
    custom_labels2 = ["A", "I", "S"]                   
    P = Array(VectorOfArray(turing_params.reaction_params))
    P2 = Array(VectorOfArray(turing_params.diffusion_constants))

    ranges = zeros(Float64, 3, 9)   #
    all_vectors = Vector{Vector{Float64}}(undef, 10)    
    for idx_j in [1;2;4;5;6;8]  #log varying parameters 
    
    for idx_i in 1:length(custom_labels)
    reaction_name = custom_labels[idx_i]
    params.reaction[reaction_name] = [P[idx_i, index]]
    end
    for idx_i in 1:length(custom_labels2)
    reaction_name = custom_labels2[idx_i]
    params.diffusion[reaction_name] = [P2[idx_i, index]]
    end
    
    
    reaction_name2 = custom_labels[idx_j]
    params.reaction[reaction_name2] = vcat(P[idx_j, index], 10 .^ range(-1.5, stop=1.5, length=70))
    @time turing_params = returnTuringParams(model, params)
    values = VectorOfArray(turing_params.reaction_params)[idx_j, :]
    all_vectors[idx_j] = values
    pmin = minimum(VectorOfArray(turing_params.reaction_params)[idx_j, :])
    pmax = maximum(VectorOfArray(turing_params.reaction_params)[idx_j, :])

    ranges[:,idx_j] = [pmin P[idx_j,index] pmax]
    end
    #for linear varying parameters (have to specify below) 
    for idx_j in [3;7]
    
    for idx_i in 1:length(custom_labels)
    reaction_name = custom_labels[idx_i]
    params.reaction[reaction_name] = [P[idx_i, index]]
    end

    for idx_i in 1:length(custom_labels2)
    reaction_name = custom_labels2[idx_i]
    params.diffusion[reaction_name] = [P2[idx_i, index]]
    end

    reaction_name2 = custom_labels[idx_j]
    params.reaction[reaction_name2] = vcat(P[idx_j, index],collect(-5:5))
    @time turing_params = returnTuringParams(model, params)
    values = VectorOfArray(turing_params.reaction_params)[idx_j, :]
    all_vectors[idx_j] = values
    pmin = minimum(VectorOfArray(turing_params.reaction_params)[idx_j, :])
    pmax = maximum(VectorOfArray(turing_params.reaction_params)[idx_j, :])

    ranges[:,idx_j] = [pmin P[idx_j,index] pmax]
    end
    # for diffusion
    for idx_j in 2:2 #check: params.diffusion 
    
    for idx_i in 1:length(custom_labels)
    reaction_name = custom_labels[idx_i]
    params.reaction[reaction_name] = [P[idx_i, index]]
    end

    for idx_i in 1:length(custom_labels2)
    reaction_name = custom_labels2[idx_i]
    params.diffusion[reaction_name] = [P2[idx_i, index]]
    end

    reaction_name2 = custom_labels2[idx_j]
    params.diffusion[reaction_name2] = vcat(P[idx_j, index],10 .^ range(-1.5, stop=1.5, length=70))
    @time turing_params = returnTuringParams(model, params)
    values = VectorOfArray(turing_params.diffusion_constants)[idx_j, :]
    all_vectors[9] = values                          #
    pmin = minimum(VectorOfArray(turing_params.diffusion_constants)[idx_j, :])
    pmax = maximum(VectorOfArray(turing_params.diffusion_constants)[idx_j, :])

    ranges[:,9] = [pmin P2[2,index] pmax]         #!!
    end
    #Plot:
    custom_labels_plot = ["k₂";"k₃";"n₁";"k₊"; "k₅"; "K₂"; "n₂";"k₆";"D_I"]   # names
    p = Plots.plot(xaxis=:log, xlims=(0.02, 40), ylims=(0,11), legend=:outerbottomright)

    for i in [1,5,6,9] # Which ones to plot
    a = all_vectors[i]
    d = ranges[:,i]
    x_midpoint =((d[1] + d[3]) / 2)
    b = size(all_vectors[i])[1]
    c = i*ones(Float64, b)
    Plots.scatter!(p, a, c, label=false,ms=3, color = "blue")
    annotate!(x_midpoint, i-0.3, custom_labels_plot[i], :color)
    Plots.scatter!(p, [d[2]], [i], label="", ms = 5, color = "red")
    Plots.title!("Robustness barcode")
    end
    display(p)

    #Plot ranges as well
    q = Plots.plot(xaxis=:log, xlims=(0.02, 40), ylims=(0,12), legend=:outerbottomright)
    for i in [1,5,6,9]   #Which ranges to plot
    a = ranges[:,i]
    x_midpoint =((a[1] + a[3]) / 2)
   
    Plots.plot!(q, a, [i,i,i], label=false,linewidth=3, color = "blue")
    annotate!(x_midpoint, i-0.3, custom_labels_plot[i], :color)
    Plots.scatter!(q, [a[2]], [i], label="", ms = 5, color = "red")
    Plots.title!("Robustness barcode")
    end
    display(q)

end

plot_robustness_barcode(turing_params, params, model, index)
#savefig("IrreversibleSequestrationRobustnessBarcode.svg")