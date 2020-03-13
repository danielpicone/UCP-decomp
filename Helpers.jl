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
        regions = subproblem.regions
        generation_vars = Dict()
        for (gen_name, gen_info) in subproblem.inputs["generators"]
            generation_vars[gen_name] = Dict()
            generation_vars[gen_name]["generation"] = Dict()
            generation_vars[gen_name]["on"] = Dict()
            for t=subproblem.start:subproblem.finish
                generation_vars[gen_name]["generation"][t] = value.(generator_output[gen_name,t])
                generation_vars[gen_name]["on"][t] = value.(generator_on[gen_name,t])
            end
        end
        df = DataFrames.DataFrame(region = String[], generator_name = String[], cost = Float64[], max_capacity = Float64[], interval = Int64[], generation = Float64[], demand = Float64[])
        for (gen_name, gen_info) in subproblem.inputs["generators"]
            r = gen_info["region"]
            for t=subproblem.start:subproblem.finish
                push!(df, [r gen_name gen_info["cost"] gen_info["capacity"] t gen_info["capacity"]*generation_vars[gen_name]["generation"][t] subproblem.demand[r][t]])
            end
        end
    return df
    end

    df = create_variables_dict(subproblem_soln)
    demand_df = create_demand_df(subproblem_soln)
    demand_df[:interval] = demand_df[:interval] .- 0.5

    # println(filter(row -> row[:interval] <= 19 & row[:interval] >= 17, df))
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
        mu[(k,k+1)] = Dict()
        for gen in keys(inputs["generators"])
            mu[(k,k+1)][gen] = Dict("ramp_up" => 0, "ramp_down" => 0)
        end
    end
    return mu
end


end
