const n_params = length(model.ps)
const n_species = length(model.states)
# convert reaction network to ODESystem
odesys = convert(ODESystem, model)

# build jacobian function
jac = ModelingToolkit.generate_jacobian(odesys,expression = Val{false})[1] #false denotes function is compiled, world issues fixed
J = ModelingToolkit.eval(jac)
jacobian(u,p,t) = J(u,p,t)

# build ODE function
f_gen = ModelingToolkit.generate_function(odesys,expression = Val{false})[1] #false denotes function is compiled, world issues fixed
f_oop = ModelingToolkit.eval(f_gen)
f_ode(u,p,t) = f_oop(u,p,t)

# Build optimized Jacobian and ODE functions using Symbolics.jl

@variables uₛₛ[1:n_species]
@parameters pₛₛ[1:n_params]

duₛₛ = Symbolics.simplify.(f_ode(collect(uₛₛ),collect(pₛₛ),0.0))

fₛₛ = eval(Symbolics.build_function(duₛₛ,vec(uₛₛ),pₛₛ;
            parallel=Symbolics.SerialForm(),expression = Val{false})[2]) #index [2] denotes in-place, mutating function
jacₛₛ = Symbolics.sparsejacobian(vec(duₛₛ),vec(uₛₛ))
fjacₛₛ = eval(Symbolics.build_function(jacₛₛ,vec(uₛₛ),pₛₛ,
            parallel=Symbolics.SerialForm(),expression = Val{false})[2]) #index [2] denotes in-place, mutating function






# build PDE function
function f_reflective(u,q)
    du_rxn = mapslices(u,dims=2) do x
        f_ode(x,q[1:n_params],0.0)
    end
    du_rxn + q[(n_params + 1):(n_params + n_species)]' .* (n_gridpoints^2*M * u)
end
   
# Build optimized Jacobian and ODE functions using Symbolics.jl

@variables uᵣ[1:n_gridpoints,1:n_species]
@parameters qᵣ[1:(n_params+n_species)]

duᵣ = Symbolics.simplify.(f_reflective(collect(uᵣ),collect(qᵣ)))



fᵣ = eval(Symbolics.build_function(duᵣ,vec(uᵣ),qᵣ;
            parallel=Symbolics.SerialForm(),expression = Val{false})[2]) #index [2] denotes in-place, mutating function
jacᵣ = Symbolics.sparsejacobian(vec(duᵣ),vec(uᵣ))
fjacᵣ = eval(Symbolics.build_function(jacᵣ,vec(uᵣ),qᵣ,
            parallel=Symbolics.SerialForm(),expression = Val{false})[2]) #index [2] denotes in-place, mutating function



function returnParam(ps,ics, i)
    n_params = length(ps)
    n_species = length(ics)
    p = zeros(Float64,n_params)
    u₀ = zeros(Float64,n_species)
    ind = CartesianIndices(Tuple([length.(ics); length.(ps)]))[Int(i)]
    for j in 1:n_species
        u₀[j] = ics[j][ind[j]]
    end
    for j in 1:n_params
        p[j] = ps[j][ind[j+n_species]]
    end
    return p,u₀
end

function returnD(ds, i)
    n_species = length(ds)
    D = zeros(Float64,n_species)
    ind = CartesianIndices(Tuple(length.(ds)))[Int(i)]
    for j in 1:n_species
        D[j] = ds[j][ind[j]]
    end
    return D
end

function computeStability(J)
    # compute general condition for Turing instability, from: https://doi.org/10.1103/PhysRevX.8.021071
     # Check if J contains NaN or Inf values                                                 # added by Daniel
   # if any(x -> isnan(x) || isinf(x), J)                                                   #added by Daniel
   #     return 0  # Indicate that J contains NaN or Inf values                                  #added by Daniel
  #  end                                                                                          #added by Daniel
    n_species = size(J)[1]
    if maximum(real(eigvals(J))) > 0.0 
        return 0  ## steady state is unstable
    else
        for subset in combinations(1:n_species)            
            if ((-1)^length(subset)*det(J[subset,subset])) < 0.0 
                 return 1 ## J is not a P_0 matrix
            end
        end
    end
    return 0 ## J is a P_0 matrix
end


 function isTuring(J,D,q2_input)

    # determine whether any eigenvalue has a real positive part across range of wavevectors (q2_input)
    #       - if yes, then store relevant variables

    max_eig = zeros(length(q2_input))

    Threads.@threads for i in eachindex(q2_input)
        max_eig[i] = maximum(real(eigvals(J - diagm(q2_input[i]*D))))
    end

    index_max = findmax(max_eig)[2]

    if findmax(max_eig)[1] > 0 & index_max < nq 
        if max_eig[length(q2_input)] < 0
            real_max = max_eig[index_max]
            q2max = sqrt(q2_input[index_max])
            qmax = sqrt(q2max)
            eigenVals = eigen(J - diagm(q2_input[index_max]*D))
            phase = sign.(eigenVals.vectors[:,findmax(real(eigenVals.values))[2]])
            if real(phase[1]) < 0
                phase = phase * -1
            end
           # M = eigvals(J - diagm(q2max*D)) 
            non_oscillatory =  Bool(imag(M[findmax(real(M))[2]]) == 0) 
            
            return qmax, phase, real_max , non_oscillatory
        else
            return 0,0,0,0 
        end
    else
        return 0,0,0,0
    end
end



function identifyTuring(sol, ds)

    # Returns turingParameters that satisfy diffusion-driven instability for each steady state in sol, and diffusion constants in ds. 
    #       - idx_ps refers to the parameter combinations (steady states [sol] are computed across parameter combinations, and do not depent on ds)
    #       - idx_ds refers to the diffusion constant combinations
    #       - idx_turing refers to the  combinations


    idx_turing = 1
    turingParameters = Array{save_turing, 1}(undef, 0)
    for solᵢ in sol
        if SciMLBase.successful_retcode(solᵢ) && minimum(solᵢ) >= 0
            J = jacobian(Array(solᵢ),solᵢ.prob.p,0.0)
            if computeStability(J) == 1
                for idx_ds in 1:prod(length.(ds))
                    qmax, phase, real_max, non_oscillatory = isTuring(J,returnD(ds,idx_ds),q2)
                    if qmax > 0
                        push!(turingParameters,save_turing(solᵢ,solᵢ.prob.p,returnD(ds,idx_ds),solᵢ.prob.u0,phase,2*pi/qmax, real_max, non_oscillatory,idx_turing))
                        idx_turing = idx_turing + 1
                    end
                end
            end
        end
    end
    return StructArray(turingParameters)
end

function get_params(model, turing_param)
    if length(turing_param.wavelength) > 1
        error("Please input only a single parameter set, not multiple (e.g., turing_params[1] instead of turing_params)")
    else
        param = model_parameters()
        param.reaction = Dict(zip(string.(parameters(model)), turing_param.reaction_params))
        param.diffusion = Dict(zip(chop.(string.(states(model)),head=0,tail=3), turing_param.diffusion_constants))
        param.domain_size = 3*turing_param.wavelength # default value
        param.initial_condition = Dict(zip(chop.(string.(states(model)),head=0,tail=3), turing_param.steady_state_values))
        param.initial_noise = 0.01
        param.random_seed = 0
        return param
    end
end

function returnParameterSets(model, params)
    # read in parameters (ps), diffusion constants (ds), initial conditions (ics), domain sizes (ls), random seeds (seeds) and random noise (noise)
    ps = Array{Array{Float64,1},1}(undef,0)
    for p in parameters(model)
        push!(ps,get!(params.reaction,string(p),[0])) ## default is zero
    end

    ds = Array{Array{Float64,1},1}(undef,0)
    ics = Array{Array{Float64,1},1}(undef,0)
    for state in states(model)
        push!(ds,get!(params.diffusion,chop(string(state), head=0,tail=3),[0])) ## default is zero
        push!(ics,get!(params.initial_condition,chop(string(state), head=0,tail=3),[1])) ## default is one
    end
    seeds = Int64.(params.random_seed)
    noise = params.initial_noise
    ls = params.domain_size
    return ps,ds,ics, ls, seeds, noise
end

function returnSingleParameter(model, params)
    # read in parameters (ps), diffusion constants (ds), initial conditions (ics), domain sizes (ls), random seeds (seeds) and random noise (noise)
    ps = Array{Float64,1}(undef,0)
    for p in parameters(model)
        push!(ps,get!(params.reaction,string(p),0)) ## default is zero
    end

    ds = Array{Float64,1}(undef,0)
    ics = Array{Float64,1}(undef,0)
    for state in states(model)
        push!(ds,get!(params.diffusion,chop(string(state), head=0,tail=3),0)) ## default is zero
        push!(ics,get!(params.initial_condition,chop(string(state), head=0,tail=3),1)) ## default is one
    end
    seeds = Int64(params.random_seed)
    noise = params.initial_noise
    ls = params.domain_size
    return ps,ds,ics, ls, seeds, noise
end

function get_param(model, turing_params, name, type)
    output = Array{Float64,1}(undef,0)
    if type == "reaction"
        labels = []
        for p in parameters(model)
            push!(labels,string(p))
        end
        index = findall(labels .== name)
        if length(index) == 0
            error("Please be sure to enter the correct parameter name and/or parameter type (reaction or diffusion); it should match the original definition in model")
        end
        for i in eachindex(turing_params)
            push!(output,turing_params[i].reaction_params[index][1])
        end
    elseif type == "diffusion"
        labels = []
        for state in states(model)
            push!(labels,chop(string(state), head=0,tail=3))
        end
        index = findall(labels .== name)
        if length(index) == 0
            error("Please be sure to enter the correct parameter name and/or parameter type (reaction or diffusion); it should match the original definition in model")
        end
        for i in eachindex(turing_params)
            push!(output,turing_params[i].diffusion_constants[index][1])
        end
    end
    return output
end

function returnTuringParams(model, params; maxiters = 1e3,alg=Rodas5(),abstol=1e-8, reltol=1e-5, tspan=1e4,ensemblealg=EnsembleThreads())   #Fills up figure better, used for figure 1 and 3
#function returnTuringParams(model, params; maxiters = 1e3,alg=Rodas5(),abstol=1e-8, reltol=1e-8, tspan=1e4,ensemblealg=EnsembleThreads())   #Daniel changed reltol from 1e-6 to 1e-8, Originally: abstol=1e-8, reltol=1e-6

    # read in parameters (ps), diffusion constants (ds), and initial conditions (ics)
   ps,ds,ics, _, _, _ = returnParameterSets(model, params)
    n_params = length(parameters(model))
    n_species = length(states(model))

 


    # Construct function to screen through parameters
    function prob_func(prob,i,repeat)
      tmp1, tmp2 = returnParam(ps,ics,i)
      @. prob.u0 = tmp2
     @. prob.p = tmp1
      return prob
    end
    
    # Define and solve EnsembleProblem to find steady state values

    p = zeros(Float64,length(ps))
    u₀ = zeros(Float64,length(ics))

    prob_fn = ODEFunction((du,u,p,t)->fₛₛ(du,vec(u),p), jac = (du,u,p,t) -> fjacₛₛ(du,vec(u),p), jac_prototype = similar(jacₛₛ,Float64))

    prob = SteadyStateProblem(prob_fn,u₀,p)
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

    sol = solve(ensemble_prob,maxiters=maxiters,DynamicSS(alg; tspan=tspan),ensemblealg,trajectories=prod([length.(ps); length.(ics)]),verbose=false,abstol=abstol, reltol=reltol, save_everystep=false)
      #Put back above trajectories=prod([length.(ps); length.(ics)])
    # Determine whether the steady state undergoes a diffusion-driven instability
    turing_params = identifyTuring(sol, ds)

    # Output simple diagnostics
    n_iter = prod([length.(ps); length.(ics); length.(ds)])
  show(stdout, MIME"text/plain"(), string("Robustness: ", round(100*length(turing_params)/n_iter),"%  (",length(turing_params),"/",n_iter,")"))

    return turing_params
end




function createIC(ic, seed, noise)
    if seed > 0
        Random.seed!(round(seed))
    end
    output = ones(n_gridpoints,length(ic))
    @. output = (ic' .* output) * abs(noise*randn() + 1)
end

function simulate(model,param;alg=KenCarp4(),reltol=1e-6,abstol=1e-8, dt = 0.1, maxiters = 1e3, save_everystep = true)  
    n_params = length(parameters(model))
    n_species = length(states(model))

    p,d,ic, l, seed, noise = returnSingleParameter(model, param)

    q_ = [p; d/(l^2)]
    u0 = createIC(ic, seed, noise)
    du = similar(u0)

 
                
    prob_fn = ODEFunction((du,u,q,t)->fᵣ(du,vec(u),q), jac = (du,u,q,t) -> fjacᵣ(du,vec(u),q), jac_prototype = similar(jacᵣ,Float64))

    ## code below is for prescribed tspan values as input argument
    # prob = ODEProblem(prob_fn,u0,tspan,q_)
    # sol = solve(prob,alg; dt=dt,reltol=reltol,abstol=abstol, maxiters=maxiters,save_everystep=save_everystep,verbose=false)
    # return sol

    prob = SteadyStateProblem(prob_fn,u0,q_)
    sol = solve(prob,DynamicSS(alg); dt=dt,reltol=reltol,abstol=abstol, maxiters=maxiters,save_everystep=save_everystep,verbose=false)

   return sol.original

end

