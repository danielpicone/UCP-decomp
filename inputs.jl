# Solving the unit commitment problem
# using CPLEX
using GLPK
using JuMP
using YAML
using DataFrames
using Gadfly
using Random


include("Helpers.jl")
include("Decomposition.jl")

# using .Decomposition
import .Decomposition
using .Helpers

# Extend optimize! for simpler syntax
import JuMP.optimize!
function optimize!(problem::Decomposition.Subproblem)
    return optimize!(problem.model)
end
function optimize!(array::Array{Decomposition.Subproblem,1})
    [optimize!(p) for p in array]
end

Random.seed!(100)

T = 24*2*2
inputs= YAML.load(open("generators-regions.yml"))
regions = Set([gen_info["region"] for (gen_name, gen_info) in inputs["generators"]])
demand = Dict()
directions = ["import", "export"]
for t=1:T
    rand_addition = rand()
    for r in regions
        if t==1
            demand[r] = Dict()
        end
        # demand[r][t] = 45+40*(rand()-0.5)*rand_addition
        demand[r][t] = 45+70*(rand()-0.5)*rand_addition
    end
end

# demand["NSW"][24] = demand["NSW"][24]-15
# demand["NSW"][25] = demand["NSW"][25]-30
# demand["NSW"][26] = demand["NSW"][26]-15

# Set up the master problem
mu = Helpers.reset_mu()

day_one = Decomposition.Subproblem(1, T, nothing, demand, inputs, regions, mu, nothing)

Decomposition.create_model!(day_one)

# Solve the model
# optimize!(day_one.model)
# Helpers.graph_subproblem(day_one)

S = 2
mu = Dict()
for s=1:S-1
    mu[(s,s+1)] = Dict()
    for gen in keys(inputs["generators"])
        mu[(s,s+1)][gen] = Dict("ramp_up" => -10000, "ramp_down" => 0)
    end
end

# day_one2 = copy(day_one)
# day_one2.mu = mu[(1, 2)]
# p1 = Decomposition.create_model!(day_one2)

# optimize!(day_one2)
# p2 = Helpers.graph_subproblem(day_one2)




# output = value.(day_one.model.obj_dict[:generator_output])
output = day_one.model.obj_dict[:generator_output]
subproblems = Decomposition.split_problem(day_one, S)
Decomposition.create_model!(subproblems)
optimize!(subproblems)


# Price the deficit subproblems
deficit_values = Dict()
for subproblem in subproblems
    s = subproblem.order
    deficit_values[s] = 0
    for trace in values(subproblem.demand)
        deficit_values[s] += sum(14700*value for value in values(trace))
    end
end


# K = 2
# rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))
#
# @variable(rmp, 1 >= lambda_deficit[s=1:S] >= 0)
# @variable(rmp, 1 >= lambda[s=1:S, 1] >= 0)
# @objective(rmp, Min, sum(deficit_values[s] * lambda_deficit[s] for s=1:S) +
#                         sum(objective_value(subproblems[s].model)*lambda[s, 1] for s=1:S))
#
# @constraint(rmp, convexity_constraint[s=1:S], sum(lambda_deficit[s]) + sum(lambda[s, k] for k=1:K) == 1)
# @constraint(rmp, ramp_down[s=1, gen in keys(subproblems[s].inputs["generators"])],
#             subproblems[s].inputs["generators"][gen]["capacity"]*sum(lambda[s, k] for k=1:K) -
#             subproblems[s].inputs["generators"][gen]["capacity"]*sum(lambda[s+1, k] for k=1:K) <=
#             subproblems[s].inputs["generators"][gen]["ramp"])
#
# @constraint(rmp, ramp_up[s=1, gen in keys(subproblems[s].inputs["generators"])],
#             subproblems[s].inputs["generators"][gen]["capacity"]*sum(lambda[s, k] for k=1:K) -
#             subproblems[s].inputs["generators"][gen]["capacity"]*sum(lambda[s+1, k] for k=1:K) >=
#             -subproblems[s].inputs["generators"][gen]["ramp"])
#
# optimize!(rmp)
# value.(lambda)


## Create the subproblems
# subproblems = []
# subproblem_solution = []
# subproblem_objective_values = []
# # convexity_duals = zeros(S)
# S = 1

# voll_cost = Array{Float64}(undef, S)
# for k=1:S
#     voll_cost[k] = sum.(v[:,k]*generators[region]["generators"]["dsm_voll"]["cost"] for (region, v) in demand) |> sum
# end
#
