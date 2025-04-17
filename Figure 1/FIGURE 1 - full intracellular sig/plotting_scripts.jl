function plot_pattern(model, sol)
    pattern = last(sol)
    x = range(0,1, length = size(pattern,1))
    labels = []
    for state in model.states
        push!(labels,chop(string(state), head=0,tail=3))
    end
    labels = reshape(labels,1,length(labels))

    pattern = pattern ./ maximum(pattern,dims=1)


    plot(x,pattern, label = labels, ticks=:none,leg = Symbol(:outer,:right), linewidth=2)


end

function plot_dynamics(model, sol; time_factor = 0.5)

    x = range(0,1, length = size(last(sol),1))    
    labels = []
    for state in model.states
        push!(labels,chop(string(state), head=0,tail=3))
    end
    labels = reshape(labels,1,length(labels))
    tfinal = sol.t[Int(ceil(time_factor*length(sol.t)))]
    tspan = range(0,tfinal,100)
    normalization_factor = maximum(last(sol),dims=1)
    for t in tspan
        normalization_factor = max(vec(normalization_factor),vec(maximum(sol(t),dims=1)))
    end
    normalization_factor = reshape(normalization_factor,1,length(normalization_factor))

    @gif for t in tspan
        pattern = sol(t)

        pattern = pattern ./ normalization_factor

        plot(x,pattern, label = labels, ticks=:none,leg = Symbol(:outer,:right), linewidth=2)
    end fps=20

end 