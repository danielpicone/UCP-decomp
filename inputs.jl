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

S = 8
T = 20*S
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
        demand[r][t] = 45+70*(rand()-0.5)*rand_addition
    end
end

new_inputs = Master.Inputs(inputs, demand)

# Set up the master problem
mu = Helpers.reset_mu()

day_one = Decomposition.Subproblem(1, T, nothing, 1, demand, inputs, regions, mu, nothing)

Decomposition.create_model!(day_one; verbose = true)
Decomposition.add_coupling_constraints!(day_one)
# optimize!(day_one.model)


subproblem_solutions = Dict{Int64, Dict{Int64, Any}}()
for s in 1:S
    subproblem_solutions[s] = Dict()
end


subproblems = Decomposition.split_problem(day_one, S)
for k=1:25
    if k == 1
        global convexity_dual = 1e10*ones(S)
    end
    Decomposition.create_model!(subproblems; verbose = true)
    optimize!(subproblems)
    convergence_value = sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S)

    println("Iteration $k: $convergence_value")
    if convergence_value > -0.001
        println("--------- Restricted master problem has reached convergence --------")
        break
    end

    for (s, sub) in enumerate(subproblems)
        # if true
        if (objective_value(subproblems[s].model) - convexity_dual[s]) < -0.001
            subproblem_solutions[s][k] = Dict("vars" => deepcopy(sub.model.obj_dict),
                                              "objective_value" => Decomposition.find_objective_value(sub))
        end
    end

    global master = Master.Master_problem(subproblems, subproblem_solutions)
    Master.create_rmp!(master)
    optimize!(master.model)
    convexity_dual = dual.(master.model.obj_dict[:convexity_constraint])

    mu = dual(master)

    # println(mu)

    for sub in subproblems
        sub.mu = mu
    end

end

Master.create_solution!(master);

mip_gap = Master.find_objective_value(master)/objective_value(master.model) - 1
