module Helpers

using JuMP
using GLPK
using DataFrames
using Gadfly

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
        storage_in = subproblem.model.obj_dict[:storage_in]
        storage_out = subproblem.model.obj_dict[:storage_out]
        regions = subproblem.regions
        generation_vars = Dict()
        for (gen_name, gen_info) in subproblem.inputs.generators
            generation_vars[gen_name] = Dict()
            generation_vars[gen_name]["generation"] = Dict()
            generation_vars[gen_name]["on"] = Dict()
            for t=subproblem.start:subproblem.finish
                generation_vars[gen_name]["generation"][t] = value.(generator_output[gen_name,t])
                generation_vars[gen_name]["on"][t] = value.(generator_on[gen_name,t])
            end
        end
        df = DataFrames.DataFrame(region = String[], name = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        for (gen_name, gen_info) in subproblem.inputs.generators
            r = gen_info["region"]
            for t=subproblem.start:subproblem.finish
                push!(df, [r gen_name gen_info["capacity"] t max(gen_info["capacity"]*generation_vars[gen_name]["generation"][t],0)])
            end
        end
        storage_vars = Dict()
        for (stg_name, stg_info) in subproblem.inputs.storage
            storage_vars[stg_name] = Dict()
            storage_vars[stg_name]["storage_in"] = Dict()
            storage_vars[stg_name]["storage_out"] = Dict()
            for t=subproblem.start:subproblem.finish
                storage_vars[stg_name]["storage_in"][t] = value.(storage_in[stg_name,t])
                storage_vars[stg_name]["storage_out"][t] = value.(storage_out[stg_name,t])
            end
        end
        stg_df = DataFrames.DataFrame(region = String[], name = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        # stg_df = DataFrames.DataFrame(region = String[], name = String[], out_or_in = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        for (stg_name, stg_info) in subproblem.inputs.storage
            r = stg_info["region"]
            for t=subproblem.start:subproblem.finish
                # push!(stg_df, [r stg_name "in" stg_info["capacity"] t -stg_info["capacity"]*storage_vars[stg_name]["storage_in"][t]])
                # push!(stg_df, [r stg_name "out" stg_info["capacity"] t stg_info["capacity"]*storage_vars[stg_name]["storage_out"][t]])
                push!(stg_df, [r string(stg_name, "_in") stg_info["capacity"] t max(stg_info["capacity"]*storage_vars[stg_name]["storage_in"][t],0)])
                push!(stg_df, [r string(stg_name, "_out") stg_info["capacity"] t max(stg_info["capacity"]*storage_vars[stg_name]["storage_out"][t], 0)])
            end
        end
        stg_df = combine(groupby(stg_df, [:region, :name, :max_capacity, :interval]), [:generation => sum])
        stg_df = rename(stg_df, [:generation_sum => :generation])
    return [df; stg_df]
    end

    df = create_variables_dict(subproblem_soln)
    demand_df = create_demand_df(subproblem_soln)
    demand_df[:interval] = demand_df[:interval] .- 0.5

    # println(filter(row -> row[:interval] <= 19 & row[:interval] >= 17, df))
    plot(df, ygroup=:region,
       Geom.subplot_grid(
       layer(sort!(demand_df, rev=true), x=:interval,y=:value, ygroup=:region,Geom.step, Theme(default_color="black")),
       layer(sort!(df, rev=true),x=:interval,y=:generation,color=:name,ygroup=:region, Geom.bar)
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
