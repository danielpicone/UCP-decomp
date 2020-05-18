# Solving the unit commitment problem
# using CPLEX
using GLPK
using JuMP
using YAML
using DataFrames
using Gadfly
using Random
using LinearAlgebra
using ProgressMeter
# using RecursiveArrayTools

S = 20
T = 24*S

inputs= YAML.load(open("generators-regions.yml"))
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
global mu = Helpers.reset_mu()

day_one = Decomposition.Subproblem(1, T, nothing, 1, demand, inputs, regions, mu, nothing)
day_one.inputs = new_inputs

Decomposition.create_model!(day_one, new_inputs; verbose = true)
# Decomposition.add_coupling_constraints!(day_one)
# optimize!(day_one.model)


subproblem_solutions = Dict{Int64, Dict{Int64, Any}}()
for s in 1:S
    subproblem_solutions[s] = Dict()
end


subproblems = Decomposition.split_problem(day_one, S)
branching_variables = Dict(i => Dict() for i in 1:S);


variable_values = Array{Array{Tuple{CartesianIndex{2},Int64},1},1}()
[push!(variable_values, []) for _ in 1:S]

root = Master.Node(0, 0, variable_values, nothing, 99999999999, false, 1, 0)

Decomposition.create_model!(subproblems, new_inputs; verbose = true)
global master = Master.Master_problem(subproblems, subproblem_solutions)

function solve_node!(node; subproblems = subproblems, master = master, mu = mu)
    # Helpers.reset_mu(mu)
    Decomposition.constrain_subproblems!(subproblems, node)
    prog = ProgressThresh(0.001, "Maximising:")
    # prog = Progress(25, 1)
    for k=1:25
        if k == 1
            global convexity_dual = 1e10*ones(S)
        end
        Decomposition.create_model!(subproblems, new_inputs; verbose = true)
        optimize!(subproblems)
        convergence_value = sum(objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S)
        convergence_value_array = [objective_value(subproblems[s].model) - convexity_dual[s] for s in 1:S]

        # println("Iteration $k: $convergence_value")
        # ProgressMeter.update!(prog, -convergence_value)
        print(lpad(string(round(convergence_value, digits = 5)), 25), "\u1b[25D")
        # next!(prog)
        if convergence_value > -0.001
            root_id = node.id
            # println("--------- Restricted master problem for node $root_id has reached convergence --------")
            break
        end

        Decomposition.add_column!(subproblems, subproblem_solutions, master, branching_variables)

        master = Master.Master_problem(subproblems, subproblem_solutions)
        Master.create_rmp!(master, new_inputs)
        Master.constrain_master!(master, node, branching_variables)
        optimize!(master.model)
        # println(objective_function(master.model))
        # println(value.(master.model[:λ]))
        convexity_dual = dual.(master.model.obj_dict[:convexity_constraint])

        mu = dual(master)
        # duals = [[value for (_, value) in x] for (_, x) in mu]
        # duals = vcat(duals...)
        # duals = vcat(duals...)
        # println(norm(duals))

        for sub in subproblems
            sub.mu = mu
        end
    end
    Decomposition.unconstrain_subproblems!(subproblems, node)
    Master.create_solution!(master)
    node.solution = master.solution
    node.objective_value = objective_value(master.model)
end

function solve_next!(; nodes_left = nodes_left, master = master,
    branching_variables = branching_variables, nodes_solved = nodes_solved,
    node_counter = node_counter, integer_found = tree_meta.integer_found, tree_meta = tree_meta)
    this_node = nodes_left[end]
    ## Fathom bad nodes
    if (this_node.objective_value > tree_meta.best_integer) & integer_found
        println("\nNode fathomed")
        pop!(nodes_left);
        return node_counter, this_node
    end
    solve_node!(this_node)
    pop!(nodes_left);
    Master.check_integer!(this_node)
    if this_node.is_integer
        tree_meta.integer_found = true
        println("\nInteger solution found")
    end
    push!(nodes_solved, this_node)
    heuristic_value = Master.find_integer_objective_value(this_node)
    if (heuristic_value != false) & (heuristic_value < tree_meta.best_integer)
        tree_meta.best_integer = heuristic_value
        tree_meta.integer_found = true
        println("Heuristic found an integer solution with value $heuristic_value")
    end
    append!(nodes_left, Master.branch(master, this_node, branching_variables, node_counter))
    node_counter += 2
    println()

    return node_counter, this_node
end
global node_counter = 0
nodes_left = []
push!(nodes_left, root)
nodes_solved = []

tree_meta = Master.Tree_meta(1, 1, T*S, 100.0, typemax(Int64), root.objective_value, false)

print("Nodes | N Iinf | Nodes remaining | MIP gap | Lower bound | Best integer")
stopping_condition = 2
while stopping_condition != 0
    Helpers.reset_mu()
    # println()
    global node_counter, this_node = solve_next!()
    tree_meta.n_nodes_left = length(nodes_left)
    tree_meta.n_iinf = nodes_left[end].n_iinf
    if stopping_condition == 2
    # if this_node.objective_value < best_incumbent
        # best_incumbent = this_node.objective_value
        tree_meta.lower_bound = this_node.objective_value
        global stopping_condition = 1
    end
    if length(nodes_left) > 0
        tree_meta.lower_bound = minimum([node.objective_value for node in nodes_left])
    end
    ## Update best integer
    if (this_node.objective_value < tree_meta.best_integer) & (this_node.is_integer)
        tree_meta.best_integer = this_node.objective_value
        tree_meta.integer_found = true
        # break
    end

    # tree_meta.mip_gap = round(min(1, tree_meta.best_integer / tree_meta.lower_bound - 1), digits = 4)*100
    tree_meta.mip_gap = min(1, tree_meta.best_integer / tree_meta.lower_bound - 1)
    if tree_meta.integer_found
        Master.sort_breadth_first!(nodes_left)
        # println("Now doing breadth first")
    else
        # Master.sort_depth_first!(nodes_left)
        Master.sort_most_integer!(nodes_left)
        # println("Doing depth first")
    end
    # Master.sort_most_integer!(nodes_left)
    if tree_meta.mip_gap < 0.01
        break
    end
    stopping_condition = length(nodes_left)
    if stopping_condition == 0
        break
    end
    tree_meta.iteration_counter += 1
    Master.print_node(tree_meta)
end

# println("The correct objective value is:               ", Decomposition.find_objective_value(day_one))
println("Branch and price found an objective value of: ", tree_meta.best_integer)
