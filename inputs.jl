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

T = 24*2
generators = YAML.load(open("generators-regions.yml"))
regions = keys(generators) |> collect
demand = Dict()
directions = ["import", "export"]
for t=1:T
    rand_addition = rand()
    for r in regions
        if t==1
            demand[r] = Dict()
        end
        demand[r][t] = 45+40*(rand()-0.5)*rand_addition
    end
end

# demand["NSW"][T+1] = demand["NSW"][T]-30

# Set up the master problem
mu = Helpers.reset_mu()

day_one = Decomposition.Subproblem(1, T, demand, generators, mu, nothing)

Decomposition.create_model!(day_one)

# Solve the model
optimize!(day_one.model)
Helpers.graph_subproblem(day_one)

mu = Dict()
for k=1:K-1
    mu[(k,k+1)] = Dict()
    for r in regions
        mu[(k,k+1)][r] = Dict()
        for gen in keys(generators[r]["generators"])
            mu[(k,k+1)][r][gen] = Dict("ramp_up" => -10000, "ramp_down" => 0)
        end
    end
end

day_one2 = copy(day_one)
day_one2.mu = mu[(1, 2)]
Decomposition.create_model!(day_one2)

optimize!(day_one2)
Helpers.graph_subproblem(day_one2)




# output = value.(day_one.model.obj_dict[:generator_output])
output = day_one.model.obj_dict[:generator_output]
K = 2
subproblems = Decomposition.split_problem(day_one, K)
Decomposition.create_model!(subproblems)
# for i in 1:K
#     println(i)
#     Decomposition.create_model!(subproblems[i])
#     optimize!(subproblems[i])
# end


# for i in output
#     println(i)
# end
# day_one = Decomposition.Subproblem_solution(day_one, optimize!(day_one.model))

## Create the subproblems
# subproblems = []
# subproblem_solution = []
# subproblem_objective_values = []
# # convexity_duals = zeros(K)
# S = 1

# voll_cost = Array{Float64}(undef, K)
# for k=1:K
#     voll_cost[k] = sum.(v[:,k]*generators[region]["generators"]["dsm_voll"]["cost"] for (region, v) in demand) |> sum
# end
#
