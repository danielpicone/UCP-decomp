module Decomposition

using JuMP
using GLPK
using DataFrames
# module for temporal Danztig-Wolfe decomposition of unit commitment problem


export Subproblem, Master, create_model!, split_problem

mutable struct Subproblem
    start::Int64
    finish::Int64
    order::Union{Int64, Nothing}
    demand::Dict
    inputs::Dict
    regions::Set
    mu::Union{Nothing, Dict}
    model::Union{Nothing, JuMP.Model}
    Subproblem(start, finish, order, demand, inputs, regions, mu, model) = start > finish ? error("Start is after finish") : new(start, finish, order,
    demand, inputs, regions, mu, model)
end

mutable struct Master
    model::JuMP.Model
    subproblem_solutions
end

Base.copy(s::Subproblem) = Subproblem(s.start,
                                      s.finish,
                                      s.order,
                                      s.demand,
                                      s.inputs,
                                      s.regions,
                                      s.mu,
                                      s.model)

function create_model!(subproblem::Subproblem; verbose = false)
    start = subproblem.start
    finish = subproblem.finish
    inputs = subproblem.inputs
    mu = subproblem.mu
    regions = subproblem.regions
    directions = ["import", "export"]

    max_cf_gens = [name for (name, gen) in inputs["generators"] if "max_cf" in keys(gen)]
    length = finish - start

    solution = Dict()
    subproblem_model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL, presolve = false))
    # Generator variables
    if verbose
        println("Adding variables")
    end
    @variable(subproblem_model, deficit[r=regions, t=start:finish] >= 0)
    @variable(subproblem_model, surplus[r=regions, t=start:finish] >= 0)
    @variable(subproblem_model, 1>=generator_output[keys(inputs["generators"]), t=start:finish]>=0)
    @variable(subproblem_model, generator_on[keys(inputs["generators"]), t=start:finish], Bin)
    @variable(subproblem_model, generator_startup[keys(inputs["generators"]), t=start:finish], Bin)

    # Interconnector variables
    @variable(subproblem_model, 1>=interconnector[keys(inputs["interconnectors"]), t=start:finish, direction = directions]>=0)

    # Objective function
    if isnothing(mu)
        @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[gen, t] +
        d["startup"]*generator_startup[gen,t] for (gen,d) in inputs["generators"],t=start:finish) +
            sum(14700*deficit[r, t] for r in regions, t=start:finish) +
            sum(1000*surplus[r, t] for r in regions, t=start:finish))
    else
        @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[gen, t] +
        d["startup"]*generator_startup[gen,t] for (gen,d) in inputs["generators"],t=start:finish) +
            sum(14700*deficit[r, t] for r in regions, t=start:finish) +
            sum(1000*surplus[r, t] for r in regions, t=start:finish) +
            sum((mu[gen]["ramp_down"]- mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, finish] for (gen, d) in inputs["generators"]) +
            sum((-mu[gen]["ramp_down"] + mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, start] for (gen, d) in inputs["generators"]))
            # sum((-mu[r][gen]["ramp_down"] + mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, start] for r in regions, (gen, d) in generators[r]["generators"]) - convexity_dual)
    end

    # Demand constraints
    if verbose
        println("Adding constraints")
    end
    @constraint(subproblem_model, demand_constraint[r=regions,t=start:finish], sum(d["capacity"]*generator_output[gen,t] for (gen,d) in inputs["generators"] if d["region"] == r)
    - sum(i["export_capacity"]*interconnector[int, t, "export"] for (int, i) in inputs["interconnectors"] if i["from"] == r)
    + sum(i["import_capacity"]*interconnector[int, t, "import"] for (int, i) in inputs["interconnectors"] if i["to"] == r)
        + deficit[r, t] - surplus[r,t] == subproblem.demand[r][t])

    # Min gen constraints
    @constraint(subproblem_model, min_gen[gen in keys(inputs["generators"]), t=start:finish], inputs["generators"][gen]["capacity"]*generator_output[gen, t] >= inputs["generators"][gen]["mingen"]*generator_on[gen, t])
    # On variable constraint
    @constraint(subproblem_model, on_constraint[gen in keys(inputs["generators"]), t=start:finish], generator_output[gen, t] <= generator_on[gen, t])
     # Ramp rate constraints
    @constraint(subproblem_model, ramp_down[gen in keys(inputs["generators"]), t=start:finish-1], inputs["generators"][gen]["capacity"]*(generator_output[gen, t] - generator_output[gen, t+1]) <= inputs["generators"][gen]["ramp"])
    @constraint(subproblem_model, ramp_up[gen in keys(inputs["generators"]), t=start:finish-1], inputs["generators"][gen]["capacity"]*(generator_output[gen, t] - generator_output[gen, t+1]) >= -inputs["generators"][gen]["ramp"])
    # Start up constraint
    @constraint(subproblem_model, start_up[gen in keys(inputs["generators"]), t=start+1:finish], generator_startup[gen,t]>=generator_on[gen,t]-generator_on[gen,t-1])
    # Maximum capacity factor constraints
    # @constraint(subproblem_model, max_cf[gen in max_cf_gens], sum(generator_output[gen, t] for t=start:finish) <=
    #             length*inputs["generators"][gen]["max_cf"])

    # Interconnector loss constraints
    @constraint(subproblem_model, interconnector_loss[int in keys(inputs["interconnectors"]), t=start:finish], (1-inputs["interconnectors"][int]["loss"])*interconnector[int, t, "export"] == interconnector[int, t, "import"])
    subproblem.model = subproblem_model
    return true
end

function create_model!(array::Array{Subproblem, 1})
    [create_model!(p) for p in array]
end

function split_problem(problem, K)
    finish = problem.finish
    length = Int64(finish/K)
    subproblems = Array{Subproblem}(undef, K)
    for k in 1:K
        subproblem_demand = Dict()
        subproblem_start = length*(k-1)
        for (r, values) in problem.demand
            subproblem_demand[r] = Dict()
            for t in (subproblem_start + problem.start):(subproblem_start + length)
                subproblem_demand[r][t] = problem.demand[r][t]
            end
        end
        println("Iteration $k: with $subproblem_demand")
        subproblems[k] = Decomposition.Subproblem(subproblem_start+1, subproblem_start + length,
        k,
        subproblem_demand,
        problem.inputs,
        problem.regions,
        problem.mu, nothing)
    end
    return subproblems
end


function update_objective(problem::Subproblem)
    start = problem.start
    finish = problem.finish
    mu = problem.mu
    @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[gen, t] +
    d["startup"]*generator_startup[gen,t] for (gen,d) in inputs["generators"],t=start:finish) +
        sum(14700*deficit[r, t] for r in regions, t=start:finish) +
        sum(1000*surplus[r, t] for r in regions, t=start:finish) +
        sum((mu[gen]["ramp_down"]- mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, finish] for (gen, d) in inputs["generators"]) +
        sum((-mu[gen]["ramp_down"] + mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, start] for (gen, d) in inputs["generators"]))
end


# function clean_solution(problem::Subproblem)
#     results = Dict()
#     for (k, v) in problem.model.obj_dict
#         if typeof(v) <: JuMP.Containers.DenseAxisArray{VariableRef, 2}
#             true
#         end
#     end
#
# end


end
