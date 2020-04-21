module Decomposition

using JuMP
using GLPK
using DataFrames
# module for temporal Danztig-Wolfe decomposition of unit commitment problem

# include("Master.jl")
using ..Master: Inputs
using ..Master: Node

export Subproblem, Master_problem, create_model!, split_problem, create_rmp!, add_coupling_constraints!, find_objective_value

mutable struct Subproblem
    start::Int64
    finish::Int64
    order::Union{Int64, Nothing}
    num_subproblems::Int64
    demand::Dict
    # inputs::Dict
    inputs
    regions::Set
    mu::Union{Nothing, Dict}
    model::Union{Nothing, JuMP.Model}
    Subproblem(start, finish, order, num_subproblems, demand, inputs, regions, mu, model) = start > finish ? error("Start is after finish") : new(start, finish, order,
    num_subproblems, demand, inputs, regions, mu, model)
end

Base.copy(s::Subproblem) = Subproblem(s.start,
                                      s.finish,
                                      s.order,
                                      s.num_subproblems,
                                      s.demand,
                                      s.inputs,
                                      s.regions,
                                      s.mu,
                                      s.model)

function create_model!(subproblem::Subproblem, inputs::Inputs; verbose = false)
    start = subproblem.start
    finish = subproblem.finish
    order = subproblem.order
    generators = inputs.generators
    interconnectors = inputs.interconnectors
    storage = inputs.storage
    intervals = start:finish
    mu = subproblem.mu
    regions = subproblem.regions
    directions = ["import", "export"]

    max_cf_gens = [name for (name, gen) in generators if :max_cf in keys(gen)]
    length = finish - start

    # subproblem_model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL, presolve = false))
    # Generator variables
    # Objective function
    if !isnothing(mu)
        generator_output = subproblem.model[:generator_output]
        generator_startup = subproblem.model[:generator_startup]
        generator_on = subproblem.model[:generator_on]
        storage_level = subproblem.model[:storage_level]
        deficit = subproblem.model[:deficit]
        surplus = subproblem.model[:surplus]
        @objective(subproblem.model, Min, sum(d["cost"] * d["capacity"] * generator_output[gen, t] +
        d["startup"]*generator_startup[gen,t] for (gen,d) in generators, t in intervals) +
            sum(14700*deficit[r, t] for r in regions, t in intervals) +
            sum(1000*surplus[r, t] for r in regions, t in intervals) +
            # sum((mu[gen]["ramp_down"]- mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, finish] for (gen, d) in inputs["generators"]) +
            # sum((-mu[gen]["ramp_down"] + mu[gen]["ramp_up"]) * d["capacity"]*generator_output[gen, start] for (gen, d) in inputs["generators"]) +
            sum(-mu[:max_cf][gen]*generator_output[gen, t] for (gen, d) in generators, t in intervals if "max_cf" in keys(d))
            + sum(mu[:master_startup][gen][order]*generator_on[gen, finish] for (gen,d) in generators if order < subproblem.num_subproblems)
            + sum(-mu[:master_startup][gen][order-1]*generator_on[gen, start] for (gen,d) in generators if order > 1)
            + sum(s["energy"]*mu[:master_storage][stg][order-1]*storage_level[stg, start] for (stg, s) in storage if order > 1)
            + sum(-s["energy"]*mu[:master_storage][stg][order]*storage_level[stg, finish+1] for (stg, s) in storage if order < subproblem.num_subproblems))
    else
        if verbose
            println("Adding variables")
            subproblem_model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL, presolve = true))
        else
            subproblem_model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF, presolve = true))
        end
        # Demand variables
        @variable(subproblem_model, deficit[r=regions, t=intervals] >= 0)
        @variable(subproblem_model, surplus[r=regions, t=intervals] >= 0)

        # Generation variables
        @variable(subproblem_model, 1>=generator_output[keys(generators), t=intervals]>=0)
        @variable(subproblem_model, generator_on[keys(generators), t=intervals], Bin)
        @variable(subproblem_model, generator_startup[keys(generators), t=intervals], Bin)

        # Interconnector variables
        @variable(subproblem_model, 1>=interconnector[keys(interconnectors), t=intervals, direction = directions]>=0)

        # Storage variables
        @variable(subproblem_model, 1>=storage_in[keys(storage), t=intervals]>=0)
        @variable(subproblem_model, 1>=storage_out[keys(storage), t=intervals]>=0)
        @variable(subproblem_model, 1>=storage_level[keys(storage), t=start:finish+1]>=0)

        @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[gen, t] +
        d["startup"]*generator_startup[gen,t] for (gen,d) in generators,t in intervals) +
            sum(14700*deficit[r, t] for r in regions, t in intervals) +
            sum(1000*surplus[r, t] for r in regions, t in intervals))

        if verbose
            println("Adding constraints")
        end
        # Demand constraints
        @constraint(subproblem_model, demand_constraint[r=regions,t=intervals], sum(d["capacity"]*generator_output[gen,t] for (gen,d) in generators if d["region"] == r)
        - sum(i["export_capacity"]*interconnector[int, t, "export"] for (int, i) in interconnectors if i["from"] == r)
        + sum(i["import_capacity"]*interconnector[int, t, "import"] for (int, i) in interconnectors if i["to"] == r)
        + sum(s["capacity"]*storage_out[stg, t] +
            - s["capacity"]*storage_in[stg, t] for (stg, s) in storage)
            + deficit[r, t] - surplus[r,t] == inputs.demand[r][t])

        # Min gen constraints
        @constraint(subproblem_model, min_gen[gen in keys(generators), t=intervals], generators[gen]["capacity"]*generator_output[gen, t] >= generators[gen]["mingen"]*generator_on[gen, t])
        # On variable constraint
        @constraint(subproblem_model, on_constraint[gen in keys(generators), t=intervals], generator_output[gen, t] <= generator_on[gen, t])
         # Ramp rate constraints
        # @constraint(subproblem_model, ramp_down[gen in keys(inputs["generators"]), t=start:finish-1], inputs["generators"][gen]["capacity"]*(generator_output[gen, t] - generator_output[gen, t+1]) <= inputs["generators"][gen]["ramp"])
        # @constraint(subproblem_model, ramp_up[gen in keys(inputs["generators"]), t=start:finish-1], inputs["generators"][gen]["capacity"]*(generator_output[gen, t] - generator_output[gen, t+1]) >= -inputs["generators"][gen]["ramp"])
        # Start up constraint
        @constraint(subproblem_model, start_up[gen in keys(generators), t=start+1:finish], generator_startup[gen,t]>=generator_on[gen,t]-generator_on[gen,t-1])

        # Interconnector loss constraints
        @constraint(subproblem_model, interconnector_loss[int in keys(interconnectors), t=start:finish], (1-interconnectors[int]["loss"])*interconnector[int, t, "export"] == interconnector[int, t, "import"])

        # Storage constraint
        @constraint(subproblem_model, storage_energy[stg in keys(storage), t in start:finish],
            storage[stg]["energy"]*storage_level[stg, t] - storage[stg]["capacity"]*storage_out[stg, t] + storage[stg]["efficiency"]*storage[stg]["capacity"]*storage_in[stg, t] == storage[stg]["energy"]*storage_level[stg, t+1])

        if (start == minimum(inputs.intervals)) & (finish == maximum(inputs.intervals))
            @constraint(subproblem_model, storage_bookends[stg in keys(storage), t in [start, finish+1]], storage_level[stg, t] == 0.5)
        elseif start == minimum(inputs.intervals)
            @constraint(subproblem_model, storage_bookends[stg in keys(storage), t in start], storage_level[stg, t] == 0.5)
        elseif finish == maximum(inputs.intervals)
            @constraint(subproblem_model, storage_bookends[stg in keys(storage), t in finish+1], storage_level[stg, t] == 0.5)
        end
        # @constraint(subproblem_model, storage_bookends[stg in keys(storage), t in [start, finish+1]], storage_level[stg, t] == 0.5)
        # @constraint(subproblem_model, storage_bookends[stg in keys(storage), t in [start, finish+1]], storage_level[stg, t] == 0.5)
        subproblem.model = subproblem_model
    end
    return true
end

function create_model!(array::Array{Subproblem, 1}, inputs::Inputs; verbose = false)
    [create_model!(p, inputs; verbose = false) for p in array]
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
        subproblems[k] = Decomposition.Subproblem(subproblem_start+1, subproblem_start + length,
        k,
        K,
        subproblem_demand,
        problem.inputs,
        problem.regions,
        problem.mu, nothing)
    end
    return subproblems
end

function find_objective_value(subproblem::Subproblem)
    inputs = subproblem.inputs
    start = subproblem.start
    finish = subproblem.finish
    generators = inputs.generators
    solution = subproblem.model.obj_dict
    obj_value = 0
    obj_value += sum(gen_info["cost"]*gen_info["capacity"]*value.(solution[:generator_output])[gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += sum(gen_info["startup"]*value.(solution[:generator_startup])[gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += 14700*sum(value.(solution[:deficit]))
    obj_value += 1000*sum(value.(solution[:surplus]))
    return obj_value
end

function add_coupling_constraints!(problem:: Subproblem)
    start = problem.start
    finish = problem.finish
    @constraint(problem.model, max_cf[gen in [gen for (gen, info) in problem.inputs.generators if "max_cf" in keys(info)]],
    sum(problem.model.obj_dict[:generator_output][gen, :]) <= finish*problem.inputs.generators[gen]["max_cf"])
    return problem.model
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


function set_branching_variable(problem::Subproblem, fix_value, ind::CartesianIndex)
    unset_binary(problem.model[:generator_on][ind])
    fix(problem.model[:generator_on][ind], fix_value)
    return problem
end

function set_branching_variable(subproblems, node::Node)
    for (s, ind, fix_value) in node.variable_values
        subproblem = subproblems[s]
        unset_binary(subproblem.model[:generator_on][ind])
        fix(subproblem.model[:generator_on][ind], fix_value)
    end
    return subproblem
end

function solve_node(subproblems::Array{Subproblem,1}, node::Node)
    for (s, fix_values) in enumerate(node.variable_values)
        subproblem = subproblems[s]
        for ((subproblem_id, ind), fix_value) in fix_values
            unset_binary(subproblem.model[:generator_on][ind])
            fix(subproblem.model[:generator_on][ind], fix_value)
        end
        optimize!(subproblem)
        for ((subproblem_id, ind), fix_value) in fix_values
            set_binary(subproblem.model[:generator_on][ind])
            unfix(subproblem.model[:generator_on][ind])
        end
    end
end

end
