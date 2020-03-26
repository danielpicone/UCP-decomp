module Master

using JuMP
using GLPK
using DataFrames
# module for temporal Danztig-Wolfe decomposition of unit commitment problem

import JuMP.dual

export Master_problem, create_rmp!, dual

mutable struct Master_problem
    inputs
    model::Union{JuMP.Model, Nothing}
    subproblems
    subproblem_solutions
    deficit_values
    S
    start
    finish
end

function Master_problem(subproblems, subproblem_solutions)
    inputs = subproblems[1].inputs
    deficit_values = Dict()
    for subproblem in subproblems
        s = subproblem.order
        deficit_values[s] = 0
        for trace in values(subproblem.demand)
            deficit_values[s] += sum(14700*value for value in values(trace))
        end
    end
    S = maximum(keys(subproblem_solutions))
    start = subproblems[1].start
    finish = subproblems[S].finish
    return Master_problem(inputs, nothing, subproblems, subproblem_solutions, deficit_values, S, start, finish)
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
    # @variable(rmp, λ_startup[s in [(i, i+1) for i in 1:S-1], gen in keys(generators)] >= 0)
    # println([subproblem_solutions[s][k]["objective_value"] for s=1:S, k in keys(subproblem_solutions[s])])
    # for s=1:S, k in keys(subproblem_solutions[s])
    #     println(subproblem_solutions[s][k]["objective_value"])
    # end
    # @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
    #                         sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])) +
    #                         sum(gen_info["startup"]*λ_startup[(s, s+1), gen] for s=1:S-1, (gen, gen_info) in generators))
    @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
                            sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])))

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
    # @constraint(rmp, master_startup[gen in keys(generators), s=1:S-1], - sum(value.(subproblem_solutions[s][k]["vars"][:generator_startup][gen, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    # + sum(value.(subproblem_solutions[s+1][k]["vars"][:generator_startup][gen, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1]))
    # <=
    # λ_startup[(s, s+1), gen])
    master.model = rmp
    return true
end

function get_solution(master::Master_problem)
    subproblem_solutions = master.subproblem_solutions
    λ_values = value.(master.model.obj_dict[:λ]).data
    generator_output = Dict()
    for (s, sub) in subproblem_solutions
        generator_output[s] = 0*value.(subproblem_solutions[s][1]["vars"][:generator_output])
        for k in keys(sub)
            generator_output[s].data .+= λ_values[(s,k)]*value.(subproblem_solutions[s][k]["vars"][:generator_output]).data
        end
    end
    return generator_output
end

function dual(master::Master_problem)
    constraints = [name for (name, value) in master.model.obj_dict if typeof(value) <: JuMP.Containers.DenseAxisArray]
    # constraints = [:max_cf, :master_startup]
    mu = Dict()
    for con in constraints
        mu[con] = Dict()
        for gen in master.model.obj_dict[con].axes[1]
            mu[con][gen] = dual.(master.model.obj_dict[con])[gen]
        end
    end
    return mu
end

end
