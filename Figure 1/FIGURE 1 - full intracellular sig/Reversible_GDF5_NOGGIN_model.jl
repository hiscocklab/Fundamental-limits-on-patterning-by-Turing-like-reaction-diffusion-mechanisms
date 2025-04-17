## Code for Reversible GDF5/NOGGIN model - Figure 1
using Plots
using RecursiveArrayTools
using Symbolics
using Statistics

## import and load scripts
include("required_scripts_2.jl")
include("plotting_scripts.jl")

function myOwnFunction(SIG,k,K,n)
    ratio = max((SIG),0)
    hill_val = ratio^n/(ratio^n + K^n)
    return k*hill_val
end

model = @reaction_network begin
    1,                      BMP --> ∅
    k₂,                     NOG --> ∅
    myOwnFunction(S24n,1,1,n₁) - k₊*BMP*NOG + k₋*C,            ∅ --> BMP 
    myOwnFunction(S24n,k₅,K₂,n₂) - k₊*BMP*NOG + k₋*C,           ∅ --> NOG 
    k₊*BMP*NOG - k₋*C,                              ∅ --> C
   BMP-0.01*R,                                               ∅ --> R
  730.0-kphos*R*S2c-kin2*S2c+kex2*S2n-0.1*S2c, ∅ --> S2c   #Constant is 73.0/Mfactor*degradation_rate (deg rate of all intracellular species set to be 0.1)
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

params.reaction["k₂"] = screen_values(min = 0.1, max = 10, mode = "log", number = 2) 
params.reaction["k₅"] = screen_values(min = 0.1, max = 10, mode = "log", number = 2)   
params.reaction["k₊"] = screen_values(min = 1, max = 10, mode = "log", number = 3)
params.reaction["k₋"] = screen_values(min = .1, max = 10, mode = "log", number = 3)  
params.reaction["K₂"] = screen_values(min = .1, max = 5, mode = "log", number = 2) 
params.reaction["n₁"] = screen_values(min = -3, max = 0, mode = "linear", number = 8)                                      
params.reaction["n₂"] = screen_values(min = -3, max = -1, mode = "linear", number = 8)    

# other parameters taken from literature, scaled proportionally to units
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

#diffusion parameters
params.diffusion["BMP"] = [1.0] 
params.diffusion["NOG"] =  screen_values(min = 0.03, max = 700, mode = "log", number = 9) 
params.diffusion["C"] =  screen_values(min = 0.03, max = 700, mode = "log", number = 9)

## Compile functions and solvers
include("simulate_scripts_2.jl")

## Turing screen
@time turing_params = returnTuringParams(model, params);  

## Plotting
index = 72 #This turing parameter set will be analysed 
param1 = get_params(model, turing_params[index]) 
param1.domain_size = 80
sol = simulate(model,param1)


## Plot solution
plot_pattern(model,sol)


## Plot just the species of interest
last(sol)
BMP = last(sol)[:,1] 
NOG =  last(sol)[:,2]  
S24n = last(sol)[:,13]

pattern= hcat(BMP, NOG, S24n)
x = range(0,1, length = size(pattern,1))
pattern = pattern ./ maximum(pattern,dims=1)

# Plot each line separately with its own label
plot(x, pattern[:,1], label = "BMP", linewidth=3, ticks=:none)
plot!(x, pattern[:,2], label = "NOG", linewidth=2)
plot!(x, pattern[:,3], label = "S24n", linewidth=2, leg=:outerright)

# Plot mRNA values
mRNA_BMP = 1 .* S24n .^n₁[index] ./(S24n .^n₁[index] .+ 1)
mRNA_BMP = mRNA_BMP ./ maximum(mRNA_BMP,dims=1)
mRNA_NOG =k₅[index] .* S24n .^ n₂[index] ./ (S24n .^ n₁[index] .+ K₂[index] .^n₂[index])
mRNA_NOG = mRNA_NOG ./ maximum(mRNA_NOG,dims=1)
plot!(x,mRNA_BMP,label = "mRNA_BMP", ticks=:none, linewidth=2)
plot!(x,mRNA_NOG,label = "mRNA_NOG", ticks=:none,linewidth=3)
