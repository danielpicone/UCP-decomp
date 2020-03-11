module Decomposition

using JuMP
using GLPK
using DataFrames
# module for temporal Danztig-Wolfe decomposition of unit commitment problem


export Subproblem, create_model!, split_problem

mutable struct Subproblem
    start::Int64
    finish::Int64
    demand::Dict
    inputs::Dict
    mu::Union{Nothing, Dict}
    model
    Subproblem(start, finish, demand, inputs, mu, model) = start > finish ? error("Start is after finish") : new(start, finish, demand, inputs, mu, model)
end

Base.copy(s::Subproblem) = Subproblem(s.start,
                                      s.finish,
                                      s.demand,
                                      s.inputs,
                                      s.mu,
                                      s.model)

function create_model!(subproblem::Subproblem; verbose = false)
    start = subproblem.start
    finish = subproblem.finish
    generators = subproblem.inputs
    mu = subproblem.mu
    regions = keys(generators) |> collect
    directions = ["import", "export"]
    solution = Dict()
    subproblem_model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))
    # Generator variables
    if verbose
        println("Adding variables")
    end
    @variable(subproblem_model, deficit[r=regions, t=start:finish] >= 0)
    @variable(subproblem_model, 0 >= surplus[r=regions, t=start:finish] >= 0)
    @variable(subproblem_model, 1>=generator_output[r=regions, keys(generators[r]["generators"]), t=start:finish]>=0)
    @variable(subproblem_model, generator_on[r=regions, keys(generators[r]["generators"]), t=start:finish], Bin)
    @variable(subproblem_model, generator_startup[r=regions, keys(generators[r]["generators"]), t=start:finish], Bin)

    # Interconnector variables
    @variable(subproblem_model, 1>=interconnector[r=regions, keys(generators[r]["interconnectors"]), t=start:finish, direction = directions]>=0)

    # Objective function
    if isnothing(mu)
        @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[r, gen, t] + d["startup"]*generator_startup[r, gen,t] for r in regions, (gen,d) in generators[r]["generators"],t=start:finish) +
            sum(14700*deficit[r, t] for r in regions, t=start:finish) +
            sum(1000*surplus[r, t] for r in regions, t=start:finish))
    else
        @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[r, gen, t] + d["startup"]*generator_startup[r, gen,t] for r in regions, (gen,d) in generators[r]["generators"],t=start:finish) +
            sum(14700*deficit[r, t] for r in regions, t=start:finish) +
            sum(1000*surplus[r, t] for r in regions, t=start:finish) +
            sum((mu[r][gen]["ramp_down"]- mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, finish] for r in regions, (gen, d) in generators[r]["generators"]) +
            sum((-mu[r][gen]["ramp_down"] + mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, start] for r in regions, (gen, d) in generators[r]["generators"]))
            # sum((-mu[r][gen]["ramp_down"] + mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, start] for r in regions, (gen, d) in generators[r]["generators"]) - convexity_dual)
    end

    # Demand constraints
    if verbose
        println("Adding constraints")
    end
    @constraint(subproblem_model, demand_constraint[r=regions,t=start:finish], sum(d["capacity"]*generator_output[r,gen,t] for (gen,d) in generators[r]["generators"]) -
        sum(i["capacity"]*interconnector[r, int, t, "export"] for (int, i) in generators[r]["interconnectors"]) + sum(i["capacity"]*interconnector[r, int, t, "import"] for (int, i) in generators[r]["interconnectors"])
        + deficit[r, t] - surplus[r,t] == subproblem.demand[r][t])

    # Min gen constraints
    @constraint(subproblem_model, min_gen[r=regions, gen in keys(generators[r]["generators"]), t=start:finish], generators[r]["generators"][gen]["capacity"]*generator_output[r, gen, t] >= generators[r]["generators"][gen]["mingen"]*generator_on[r, gen, t])
    # On variable constraint
    @constraint(subproblem_model, on_constraint[r=regions, gen in keys(generators[r]["generators"]), t=start:finish], generator_output[r, gen, t] <= generator_on[r, gen, t])
     # Ramp rate constraints
    @constraint(subproblem_model, ramp_down[r=regions, gen in keys(generators[r]["generators"]), t=start:finish-1], generators[r]["generators"][gen]["capacity"]*(generator_output[r, gen, t] - generator_output[r, gen, t+1]) <= generators[r]["generators"][gen]["ramp"])
    @constraint(subproblem_model, ramp_up[r=regions, gen in keys(generators[r]["generators"]), t=start:finish-1], generators[r]["generators"][gen]["capacity"]*(generator_output[r, gen, t] - generator_output[r, gen, t+1]) >= -generators[r]["generators"][gen]["ramp"])
    # Start up constraint
    @constraint(subproblem_model, start_up[r=regions, gen in keys(generators[r]["generators"]), t=start+1:finish], generator_startup[r, gen,t]>=generator_on[r, gen,t]-generator_on[r, gen,t-1])

    # Interconnector loss constraints
    @constraint(subproblem_model, interconnector_loss[r=regions, int in keys(generators[r]["interconnectors"]), t=start:finish], (1-generators[r]["interconnectors"][int]["loss"])*interconnector[r, int, t, "export"] == interconnector[generators[r]["interconnectors"][int]["to"], int, t, "import"])
    # println(generator_output)
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
    for i in 1:K
        subproblem_demand = Dict()
        subproblem_start = length*(i-1)
        for (r, values) in problem.demand
            subproblem_demand[r] = Dict()
            for t in (subproblem_start + problem.start):(subproblem_start + length)
                subproblem_demand[r][t] = problem.demand[r][t]
            end
        end
        println("Iteration $i: with $subproblem_demand")
        subproblems[i] = Decomposition.Subproblem(subproblem_start+1, subproblem_start + length, subproblem_demand, problem.inputs, problem.mu, nothing)
    end
    return subproblems
end


function update_objective(problem::Subproblem)
    start = problem.start
    finish = problem.finish
    mu = problem.mu
    @objective(subproblem_model, Min, sum(d["cost"] * d["capacity"] * generator_output[r, gen, t] + d["startup"]*generator_startup[r, gen,t] for r in regions, (gen,d) in generators[r]["generators"],t=start:finish) +
        sum(14700*deficit[r, t] for r in regions, t=start:finish) +
        sum(1000*surplus[r, t] for r in regions, t=start:finish) +
        sum((mu[r][gen]["ramp_down"]- mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, finish] for r in regions, (gen, d) in generators[r]["generators"]) +
        sum((-mu[r][gen]["ramp_down"] + mu[r][gen]["ramp_up"]) * d["capacity"]*generator_output[r, gen, start] for r in regions, (gen, d) in generators[r]["generators"]))
end


end
