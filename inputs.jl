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
include("Master.jl")
#
using .Decomposition
# import .Decomposition
using .Helpers
using .Master

# Extend optimize! for simpler syntax
import JuMP.optimize!
function optimize!(problem::Decomposition.Subproblem)
    return optimize!(problem.model)
end
function optimize!(array::Array{Decomposition.Subproblem,1})
    [optimize!(p) for p in array]
end

Random.seed!(100)

S = 6
T = 24*2*S
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


subproblem_solutions = Dict{Int64, Dict{Int64, Any}}()
for s in 1:S
    subproblem_solutions[s] = Dict()
end


subproblems = Decomposition.split_problem(day_one, S)
for k=1:10
    if k == 1
        global convexity_dual = 1e10*ones(S)
    end
    Decomposition.create_model!(subproblems; verbose = true)
    optimize!(subproblems)
    # [println(objective_value(s.model)) for s in subproblems]
    # println("This is the convexity_dual $convexity_dual")
    convergence_value = sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S)
    println("Iteration $k: $convergence_value")
    if sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S) > -0.001
        break
    end

    for (s, sub) in enumerate(subproblems)
        subproblem_solutions[s][k] = Dict("vars" => sub.model.obj_dict,
                                          "objective_value" => objective_value(sub.model))
    end
    master = Master.Master_problem(subproblems, subproblem_solutions)
    Master.create_rmp!(master)

    optimize!(master.model)
    # println(value.(master.model.obj_dict[:λ]))
    # println(dual.(master.model.obj_dict[:convexity_constraint]))
    convexity_dual = dual.(master.model.obj_dict[:convexity_constraint])

    mu = dual(master)

    for sub in subproblems
        sub.mu = mu
    end

end


# subproblems[1].mu = mu
# subproblems[2].mu = mu
# Decomposition.create_model!(subproblems)
# optimize!(subproblems)
# for (i, s) in enumerate(subproblems)
#     subproblem_solutions[i][2] = Dict("vars" => s.model.obj_dict,
#                                       "objective_value" => objective_value(s.model))
# end
#
# master = Master.Master_problem(subproblems, subproblem_solutions)
# Master.create_rmp!(master)
# optimize!(master.model)
