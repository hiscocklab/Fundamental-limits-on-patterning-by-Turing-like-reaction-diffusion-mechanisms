using Catalyst, Combinatorics, Random, StructArrays
using DifferentialEquations, LinearAlgebra, ModelingToolkit, Symbolics
using JLD2

mutable struct save_turing
    steady_state_values::Vector{Float64}
    reaction_params::Vector{Float64}
    diffusion_constants::Vector{Float64}
    initial_conditions::Vector{Float64}
    pattern_phase::Vector{Int64}
    wavelength::Float64
    max_real_eigval::Float64
    non_oscillatory::Bool
    idx_turing::Int64
end


mutable struct model_parameters
    reaction
    diffusion
    initial_condition
    initial_noise
    domain_size
    random_seed
    function model_parameters()
        return new(Dict(),Dict(),Dict(),0.01,[1.0],[1])  
    end
    
end

"""
    screen_values(;min=0,max=1,number=10, mode="linear")

Return double the number `x` plus `1`.
"""
function screen_values(;min=0,max=1,number=10, mode="linear")
    if number > 1
        if mode == "linear"
            return collect((range(min,stop=max,length=number)))
        elseif mode == "log"
            return collect(10 .^(range(log10(min),stop=log10(max),length=number)))
        elseif mode == "seed"
            return collect(range(1,stop=N,length=number))
        end
    else 
        error("Please ensure number of values is greater than 1")
    end
end

function parameterSweep(model;min=0,max=1,number=10, mode="linear", var = "param")
    sweep = sweep_vals(min=min,max=max,number=number,mode=mode)
    if var == "param"
        return Dict(zip(string.(parameters(model)), [sweep for i in 1:length(parameters(model))]))
    elseif var == "diffusion"
        return Dict(zip(string.(states(model)), [sweep for i in 1:length(states(model))]))
    end
end

const nq = 100
const q2 =  10 .^ (range(-2, stop=2, length=nq))                                               #10 .^(range(-2,stop=2,length=nq))
const n_denser = (range(0,stop=1,length=100))
const n_gridpoints = 128
const M = Array(Tridiagonal([1.0 for i in 1:n_gridpoints-1],[-2.0 for i in 1:n_gridpoints],[1.0 for i in 1:n_gridpoints-1]))
M[1,2] = 2.0
M[end,end-1] = 2.0


  