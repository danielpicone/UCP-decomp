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

T = 24*2*4
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

S = 4
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

subproblem_solutions = Dict{Int64, Dict{Int64, Any}}()
for (i, s) in enumerate(subproblems)
    subproblem_solutions[i] = Dict()
    subproblem_solutions[i][1] = s.model.obj_dict
end


K = 1
rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))

@variable(rmp, 1 >= lambda_deficit[s=1:S] >= 0)
@variable(rmp, 1 >= λ[s=1:S, 1] >= 0)
@objective(rmp, Min, sum(deficit_values[s] * lambda_deficit[s] for s=1:S) +
                        sum(objective_value(subproblems[s].model)*λ[s, 1] for s=1:S))

@constraint(rmp, convexity_constraint[s=1:S], sum(lambda_deficit[s]) + sum(λ[s, k] for k=1:K) == 1)
# Ramp down constraints
@constraint(rmp, ramp_down[s=1, gen in keys(subproblems[s].inputs["generators"])],
            value(subproblems[s].model.obj_dict[:generator_output][gen, subproblems[s].finish])*sum(λ[s, k] for k=1:K) -
            value(subproblems[s+1].model.obj_dict[:generator_output][gen, subproblems[s+1].start])*sum(λ[s+1, k] for k=1:K) <=
            subproblems[s].inputs["generators"][gen]["ramp"])

# Ramp up constraints
@constraint(rmp, ramp_up[s=1, gen in keys(subproblems[s].inputs["generators"])],
            value(subproblems[s].model.obj_dict[:generator_output][gen, subproblems[s].finish])*sum(λ[s, k] for k=1:K) -
            value(subproblems[s+1].model.obj_dict[:generator_output][gen, subproblems[s+1].start])*sum(λ[s+1, k] for k=1:K) >=
            -subproblems[s].inputs["generators"][gen]["ramp"])

# Max cf constraints
@constraint(rmp, max_cf[gen in [key for (key, value) in subproblems[1].inputs["generators"] if "max_cf" in keys(value)]],
            sum( λ[s, 1]*value(subproblems[s].model.obj_dict[:generator_output][gen, t]) for s=1:S, t=subproblems[s].start:subproblems[s].finish)
            <= subproblems[S].finish*subproblems[1].inputs["generators"][gen]["max_cf"])

optimize!(rmp)
value.(λ)


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
