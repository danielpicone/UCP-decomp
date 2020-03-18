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
# function create_rmp!(master)
    subproblem_solutions = master.subproblem_solutions
    subproblems = master.subproblems
    # println("We are now in the RMP creation function")
    S = master.S
    # println(S)
    # rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))
    rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF))
    # master.model = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))

    # println(rmp)
    @variable(rmp, 1 >= lambda_deficit[s=1:S] >= 0)
    @variable(rmp, 1 >= λ[s=1:S, k in keys(subproblem_solutions[s])] >= 0)
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
                sum( λ[s, k]*value(subproblem_solutions[s][k]["vars"][:generator_output][gen, t]) for s=1:S, k in keys(subproblem_solutions[s]), t=subproblems[s].start:subproblems[s].finish)
                <= master.finish*master.inputs["generators"][gen]["max_cf"])
    master.model = rmp
    return true
end

function dual(master::Master_problem)
    constraints = [name for (name, value) in master.model.obj_dict if typeof(value) <: JuMP.Containers.DenseAxisArray]
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
