module Helpers

using JuMP
using GLPK
using DataFrames
using Gadfly
# include("Decomposition.jl")

# using .Decomposition
# Helper functions for easy problem creation

export graph_subproblem, reset_mu

function graph_subproblem(subproblem_soln)
    function create_demand_df(subproblem)
        d = DataFrames.DataFrame(region = String[], interval = Float64[], value = Float64[])
        for (r, values) in subproblem.demand
            for t in subproblem.start:subproblem.finish
                push!(d, [r t values[t]])
            end
            push!(d, [r subproblem.finish+1 values[subproblem.finish]])
        end
        return d
    end

    function create_variables_dict(subproblem)
        generator_output = subproblem.model.obj_dict[:generator_output]
        generator_on = subproblem.model.obj_dict[:generator_on]
        regions = keys(subproblem.inputs) |> collect
        generation_vars = Dict()
        for r in regions
            generation_vars[r] = Dict()
            for gen in keys(subproblem.inputs[r]["generators"])
                generation_vars[r][gen] = Dict()
                generation_vars[r][gen]["generation"] = Dict()
                generation_vars[r][gen]["on"] = Dict()
                for t=subproblem.start:subproblem.finish
                    generation_vars[r][gen]["generation"][t] = value.(generator_output[r, gen,t])
                    generation_vars[r][gen]["on"][t] = value.(generator_on[r, gen,t])
                end
            end
        end
        df = DataFrames.DataFrame(region = String[], generator_name = String[], cost = Float64[], max_capacity = Float64[], interval = Int64[], generation = Float64[], demand = Float64[])
        for r in regions
            for (gen_name, gen_info) in subproblem.inputs[r]["generators"]
                for t=subproblem.start:subproblem.finish
                    push!(df, [r gen_name gen_info["cost"] gen_info["capacity"] t gen_info["capacity"]*generation_vars[r][gen_name]["generation"][t] subproblem.demand[r][t]])
                end
            end
        end
    return df
    end

    df = create_variables_dict(subproblem_soln)
    demand_df = create_demand_df(subproblem_soln)
    demand_df[:interval] = demand_df[:interval] .- 0.5

    plot(df, ygroup=:region,
       Geom.subplot_grid(
       layer(sort!(demand_df, rev=true), x=:interval,y=:value, ygroup=:region,Geom.step, Theme(default_color="black")),
       layer(sort!(df, rev=true),x=:interval,y=:generation,color=:generator_name,ygroup=:region, Geom.bar)
       ))
end

function reset_mu(K=nothing)
    if isnothing(K)
        return nothing
    end
    mu = Dict()
    for k=1:K-1
        mu[(k,k+1)] = Dict()
        for r in regions
            mu[(k,k+1)][r] = Dict()
            for gen in keys(generators[r]["generators"])
                mu[(k,k+1)][r][gen] = Dict("ramp_up" => 0, "ramp_down" => 0)
            end
        end
    end
    return mu
end


end
