# Solving the unit commitment problem
# using CPLEX
using GLPK
using JuMP
using YAML
using DataFrames
using Gadfly
using Random

S = 3
T = 5*S
node_id = 0

include("Helpers.jl")
include("Master.jl")
include("Decomposition.jl")

using .Master
using .Decomposition
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
day_one.inputs = new_inputs

Decomposition.create_model!(day_one, new_inputs; verbose = true)
Decomposition.add_coupling_constraints!(day_one)
optimize!(day_one.model)


subproblem_solutions = Dict{Int64, Dict{Int64, Any}}()
for s in 1:S
    subproblem_solutions[s] = Dict()
end


subproblems = Decomposition.split_problem(day_one, S)
for k=1:25
    if k == 1
        global convexity_dual = 1e10*ones(S)
    end
    Decomposition.create_model!(subproblems, new_inputs; verbose = true)
    optimize!(subproblems)
    convergence_value = sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S)

    println("Iteration $k: $convergence_value")
    if convergence_value > -0.001
        println("--------- Restricted master problem has reached convergence --------")
        break
    end

    for (s, sub) in enumerate(subproblems)
        if (objective_value(subproblems[s].model) - convexity_dual[s]) < -0.001
            subproblem_solutions[s][k] = Dict("vars" => deepcopy(sub.model.obj_dict),
                                              "objective_value" => Decomposition.find_objective_value(sub))
        end
    end

    global master = Master.Master_problem(subproblems, subproblem_solutions)
    Master.create_rmp!(master, new_inputs)
    optimize!(master.model)
    convexity_dual = dual.(master.model.obj_dict[:convexity_constraint])

    mu = dual(master)

    # println(mu)

    for sub in subproblems
        sub.mu = mu
    end

end

Master.create_solution!(master);

branching_variables = Dict();

mip_gap = Master.find_integer_objective_value(master)/objective_value(master.model) - 1

# id, ind = Master.most_fractional_index(master)
# before = value.(subproblems[id].model[:generator_on])


variable_values = Array{Array{Tuple{CartesianIndex{2},Int64},1},1}()
[push!(variable_values, []) for _ in 1:S]

root = Master.Node(0, 0, variable_values, nothing)
# # root = Master.Node(0, 0, variable_values, master.solution)
# left, right = Master.branch(master, root, branching_variables)
# node_id += 2
# Decomposition.solve_node(subproblems, left)
#
# Decomposition.solve_node(subproblems, right)

function solve_root!(root)
    for k=1:25
        if k == 1
            global convexity_dual = 1e10*ones(S)
        end
        Decomposition.create_model!(subproblems, new_inputs; verbose = true)
        optimize!(subproblems)
        convergence_value = sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S)

        println("Iteration $k: $convergence_value")
        if convergence_value > -0.001
            println("--------- Restricted master problem has reached convergence --------")
            break
        end

        for (s, sub) in enumerate(subproblems)
            if (objective_value(subproblems[s].model) - convexity_dual[s]) < -0.001
                subproblem_solutions[s][k] = Dict("vars" => deepcopy(sub.model.obj_dict),
                                                  "objective_value" => Decomposition.find_objective_value(sub))
            end
        end

        global master = Master.Master_problem(subproblems, subproblem_solutions)
        Master.create_rmp!(master, new_inputs)
        optimize!(master.model)
        convexity_dual = dual.(master.model.obj_dict[:convexity_constraint])

        mu = dual(master)

        # println(mu)

        for sub in subproblems
            sub.mu = mu
        end

    end
    Master.create_solution!(master)
    root.solution = master.solution
end

solve_root!(root)
