module Master

using JuMP
using GLPK
using DataFrames
using Gadfly
# module for temporal Danztig-Wolfe decomposition of unit commitment problem

variable_names = [:generator_output, :generator_on, :deficit, :surplus]

import JuMP.dual

export Master_problem, create_rmp!, dual, create_solution!, integer_heuristic

mutable struct Master_problem
    inputs
    regions
    demand
    model::Union{JuMP.Model, Nothing}
    subproblems
    subproblem_solutions
    deficit_values
    S
    start
    finish
    solution::Union{Dict{Symbol, JuMP.Containers.DenseAxisArray}, Nothing}
end

struct Inputs
    demand
    generators
    interconnectors
    storage
    regions
    intervals
end

function Inputs(inputs::Dict, demand::Dict)
    demand = demand
    generators = get(inputs ,"generators", false)
    interconnectors = get(inputs, "interconnectors", false)
    storage = get(inputs, "storage", false)
    regions = Set([gen_info["region"] for (gen_name, gen_info) in generators])
    start = minimum([minimum(keys(trace)) for (r, trace) in demand])
    finish = maximum([maximum(keys(trace)) for (r, trace) in demand])
    intervals = start:finish

    return Inputs(demand, generators, interconnectors, storage, regions, intervals)
end


function Base.show(io::IO, master::Master_problem)
    S = master.S
    num_vars = length(all_variables(master.model))
    print(io, "RMP with $S subproblems and $num_vars variables")
end

function Master_problem(subproblems, subproblem_solutions)
    inputs = subproblems[1].inputs
    regions = subproblems[1].regions
    deficit_values = Dict()
    demand = Dict(r => Dict() for r in regions)
    for subproblem in subproblems
        s = subproblem.order
        deficit_values[s] = 0
        # for trace in values(subproblem.demand)
        for (r, trace) in subproblem.demand
            deficit_values[s] += sum(14700*value for value in values(trace))
            # println(demand[r])
            # println(trace)
            demand[r] = merge(demand[r], trace)
        end
    end
    S = maximum(keys(subproblem_solutions))
    start = subproblems[1].start
    finish = subproblems[S].finish
    return Master_problem(inputs, regions, demand, nothing, subproblems, subproblem_solutions, deficit_values, S, start, finish, nothing)
end


function create_rmp!(master::Master_problem)
    subproblem_solutions = master.subproblem_solutions
    subproblems = master.subproblems
    inputs = master.inputs
    start = master.start
    finish = master.finish
    generators = inputs["generators"]
    S = master.S
    # rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))
    rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF))

    @variable(rmp, 1 >= lambda_deficit[s=1:S] >= 0)
    @variable(rmp, 1 >= λ[s=1:S, k in keys(subproblem_solutions[s])] >= 0)
    @variable(rmp, λ_startup[s in [(i, i+1) for i in 1:S-1], gen in keys(generators)] >= 0)
    @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
                            sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])) +
                            sum(gen_info["startup"]*λ_startup[(s, s+1), gen] for s=1:S-1, (gen, gen_info) in generators))
    # @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
    #                         sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])))

    @constraint(rmp, convexity_constraint[s=1:S], sum(lambda_deficit[s]) + sum(λ[s, k] for k in keys(subproblem_solutions[s])) == 1)
    # Ramp down constraints
    # @constraint(rmp, ramp_down[s=1, gen in keys(subproblems[s].inputs["generators"])],
    #             sum(value(subproblem_solutions[s][k]["vars"][:generator_output][gen, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s])) -
    #             sum(value(subproblem_solutions[s+1][k]["vars"][:generator_output][gen, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1])) <=
    #             master.inputs["generators"][gen]["ramp"])

    # Ramp up constraints
    # @constraint(rmp, ramp_up[s=1, gen in keys(subproblems[s].inputs["generators"])],
    #             sum(value(subproblem_solutions[s][k]["vars"][:generator_output][gen, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s])) -
    #             sum(value(subproblem_solutions[s+1][k]["vars"][:generator_output][gen, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1])) >=
    #             -master.inputs["generators"][gen]["ramp"])

    # Max cf constraints
    @constraint(rmp, max_cf[gen in [key for (key, value) in master.inputs["generators"] if "max_cf" in keys(value)]],
                # sum( λ[s, k]*sum(value(subproblem_solutions[s][k]["vars"][:generator_output][gen, t]) for s=1:S, k in keys(subproblem_solutions[s]), t=subproblems[s].start:subproblems[s].finish) for s=1:S, k in keys(subproblem_solutions[s]))
                sum( λ[s, k]*value(subproblem_solutions[s][k]["vars"][:generator_output][gen, t]) for s=1:S, k in keys(subproblem_solutions[s]), t=subproblems[s].start:subproblems[s].finish)
                <= master.finish*master.inputs["generators"][gen]["max_cf"])
    # Startup cost constraints
    @constraint(rmp, master_startup[gen in keys(generators), s=1:S-1], - sum(value.(subproblem_solutions[s][k]["vars"][:generator_on][gen, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    + sum(value.(subproblem_solutions[s+1][k]["vars"][:generator_on][gen, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1]))
    <=
    λ_startup[(s, s+1), gen])
    master.model = rmp
    return true
end

function get_solution(master::Master_problem)
    subproblem_solutions = master.subproblem_solutions
    λ_values = value.(master.model.obj_dict[:λ]).data
    variables = Dict()
    for var in variable_names
        var_dict = Dict()
        for (s, sub) in subproblem_solutions
            var_dict[s] = 0*value.(subproblem_solutions[s][1]["vars"][var])
            for k in keys(sub)
                var_dict[s].data .+= λ_values[(s,k)]*value.(subproblem_solutions[s][k]["vars"][var]).data
            end
        end
        variables[var] = var_dict
    end
    return variables
end

function dual(master::Master_problem)
    # constraints = [name for (name, value) in master.model.obj_dict if typeof(value) <: JuMP.Containers.DenseAxisArray]
    constraints = [:max_cf, :master_startup]
    mu = Dict()
    for con in constraints
        mu[con] = Dict()
        for gen in master.model.obj_dict[con].axes[1]
            if con == :max_cf
                mu[con][gen] = dual.(master.model.obj_dict[con])[gen]
            elseif con == :master_startup
                mu[con][gen] = dual.(master.model.obj_dict[con])[gen,:].data
            end
        end
    end
    return mu
end

function create_solution!(master::Master_problem)
    master_solutions = get_solution(master)
    num_subproblems = master.S
    master.solution = Dict()
    for var in variable_names
        new_data = hcat([master_solutions[var][i].data for i in 1:num_subproblems]...) .|>
                        x->round(x; digits = 10) |> x->maximum([x,0])

        gens = master_solutions[var][1].axes[1]
        intervals = master.start:master.finish
        array = JuMP.Containers.DenseAxisArray(new_data, gens, intervals)
        master.solution[var] = array
    end
    return true
end

function graph_master(master::Master_problem)
    function create_demand_df(master)
        d = DataFrames.DataFrame(region = String[], interval = Float64[], value = Float64[])
        for (r, values) in master.demand
            for t in master.start:master.finish
                push!(d, [r t values[t]])
            end
            push!(d, [r master.finish+1 values[master.finish]])
        end
        return d
    end

    function create_variables_dict(master)
        generator_output = master.solution[:generator_output]
        regions = master.regions
        generation_vars = Dict()
        for (gen_name, gen_info) in master.inputs["generators"]
            generation_vars[gen_name] = Dict()
            generation_vars[gen_name]["generation"] = Dict()
            for t=master.start:master.finish
                generation_vars[gen_name]["generation"][t] = generator_output[gen_name,t]
            end
        end
        df = DataFrames.DataFrame(region = String[], generator_name = String[], cost = Float64[], max_capacity = Float64[], interval = Int64[], generation = Float64[], demand = Float64[])
        for (gen_name, gen_info) in master.inputs["generators"]
            r = gen_info["region"]
            for t=master.start:master.finish
                push!(df, [r gen_name gen_info["cost"] gen_info["capacity"] t max(gen_info["capacity"]*generation_vars[gen_name]["generation"][t],0) master.demand[r][t]])
            end
        end
    return df
    end

    df = create_variables_dict(master)
    demand_df = create_demand_df(master)
    demand_df[:interval] = demand_df[:interval] .- 0.5

    plot(df, ygroup=:region,
       Geom.subplot_grid(
       layer(sort!(demand_df, rev=true), x=:interval,y=:value, ygroup=:region,Geom.step, Theme(default_color="black")),
       layer(sort!(df, rev=true),x=:interval,y=:generation,color=:generator_name,ygroup=:region, Geom.bar)
       ))
end

function integer_heuristic(master)
    # First make generator_on variables binary
    axes = master.solution[:generator_output].axes
    integer_solution = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:generator_output => master.solution[:generator_output])
    integer_generator_on = ceil.(master.solution[:generator_on])
    integer_solution[:generator_on] = integer_generator_on
    num_gens = size(integer_generator_on)[1]

    # Then calculate the new generator_startup values
    integer_generator_startup = ([integer_generator_on zeros(num_gens,1)] - [zeros(num_gens,1) integer_generator_on])[:, 1:end-1]
    integer_solution[:generator_startup] = max.(JuMP.Containers.DenseAxisArray(integer_generator_startup, axes...),0)

    integer_solution[:deficit] = master.solution[:deficit]
    integer_solution[:surplus] = master.solution[:surplus]

    return integer_solution
end

function find_objective_value(master::Master_problem)
    integer_solution = integer_heuristic(master)

    inputs = master.inputs
    start = master.start
    finish = master.finish
    generators = inputs["generators"]
    obj_value = 0
    obj_value += sum(gen_info["cost"]*gen_info["capacity"]*integer_solution[:generator_output][gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += sum(gen_info["startup"]*integer_solution[:generator_startup][gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += 14700*sum(integer_solution[:deficit])
    obj_value += 1000*sum(integer_solution[:surplus])
    return obj_value
end


end
