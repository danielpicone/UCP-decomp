module Master

using JuMP
using GLPK
using DataFrames
using Gadfly
# module for temporal Danztig-Wolfe decomposition of unit commitment problem

variable_names = [:generator_output, :generator_on, :generator_startup, :deficit, :surplus, :storage_in, :storage_out]
INDEX_TYPE = Tuple{Any, Any}
# INDEX_TYPE = Tuple{CartesianIndex{2},Int64}

import JuMP.dual
using Main: S, T, inputs

export Master_problem, Tree, Inputs, Node, Tree_meta,
       create_rmp!, dual, create_solution!, integer_heuristic,
       add_branching_variable!

mutable struct Master_problem
    inputs
    regions
    demand
    model::Union{JuMP.Model, Nothing}
    subproblems
    subproblem_solutions
    deficit_values
    S
    start
    finish
    solution::Union{Dict{Symbol, JuMP.Containers.DenseAxisArray}, Nothing}
end

struct Inputs
    demand
    generators
    interconnectors
    storage
    regions
    intervals
end

mutable struct Node
    id::Int32
    mother_node_id::Int32
    variable_values::Array{Array{INDEX_TYPE,1},1}
    solution::Union{Dict{Symbol, JuMP.Containers.DenseAxisArray}, Nothing}
    objective_value
    is_integer::Bool
    n_iinf::Int64
    depth::Int32
end

mutable struct Tree
    variables::Dict{Int64, Tuple}
    depth::Int64
end

mutable struct Tree_meta
    iteration_counter::Int64
    n_nodes_left::Int64
    n_iinf::Int32
    mip_gap::Float64
    best_integer
    lower_bound
    integer_found::Bool
end

function Inputs(inputs::Dict, demand::Dict)
    demand = demand
    generators = get(inputs ,"generators", false)
    interconnectors = get(inputs, "interconnectors", false)
    storage = get(inputs, "storage", false)
    regions = Set([gen_info["region"] for (gen_name, gen_info) in generators])
    start = minimum([minimum(keys(trace)) for (r, trace) in demand])
    finish = maximum([maximum(keys(trace)) for (r, trace) in demand])
    intervals = start:finish

    return Inputs(demand, generators, interconnectors, storage, regions, intervals)
end

function Node(master::Master_problem, id, mother_node::Node, variable_key, fix_value)
    new_variable_values = deepcopy(mother_node.variable_values)
    append!(new_variable_values[variable_key[1]], [(variable_key, fix_value)])
    return Node(id, mother_node.id, new_variable_values, nothing, mother_node.objective_value, false, mother_node.n_iinf-1, mother_node.depth+1)
end

function Master_problem(subproblems, subproblem_solutions)
    inputs = subproblems[1].inputs
    regions = subproblems[1].regions
    deficit_values = Dict()
    demand = Dict(r => Dict() for r in regions)
    for subproblem in subproblems
        s = subproblem.order
        deficit_values[s] = 0
        for (r, trace) in subproblem.demand
            deficit_values[s] += sum(14700*value for value in values(trace))
            demand[r] = merge(demand[r], trace)
        end
    end
    S = maximum(keys(subproblems))
    start = subproblems[1].start
    finish = subproblems[S].finish
    if isnothing(subproblem_solutions)
        subproblem_solutions = Dict(i => Dict() for i in 1:S)
    end
    return Master_problem(inputs, regions, demand, nothing, subproblems, subproblem_solutions, deficit_values, S, start, finish, nothing)
end


function create_rmp!(master::Master_problem, inputs::Inputs)
    subproblem_solutions = master.subproblem_solutions
    subproblems = master.subproblems
    start = minimum(inputs.intervals)
    finish = maximum(inputs.intervals)
    generators = inputs.generators
    storage = inputs.storage
    S = master.S
    # rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))
    rmp = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.OFF))
    # rmp = JuMP.Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.OFF))

    @variable(rmp, 1 >= lambda_deficit[s=1:S] >= 0)
    @variable(rmp, 1 >= λ[s=1:S, k in keys(subproblem_solutions[s])] >= 0)
    @variable(rmp, λ_startup[s in [(i, i+1) for i in 1:S-1], gen in keys(generators)] >= 0)
    @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
                            sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])) +
                            sum(gen_info["startup"]*λ_startup[(s, s+1), gen] for s=1:S-1, (gen, gen_info) in generators))
    # @objective(rmp, Min, sum(master.deficit_values[s] * lambda_deficit[s] for s=1:S) +
    #                         sum(subproblem_solutions[s][k]["objective_value"]*λ[s, k] for s=1:S, k in keys(subproblem_solutions[s])))

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
    @constraint(rmp, max_cf[gen in [key for (key, value) in generators if "max_cf" in keys(value)]],
                # sum( λ[s, k]*sum(value(subproblem_solutions[s][k]["vars"][:generator_output][gen, t]) for s=1:S, k in keys(subproblem_solutions[s]), t=subproblems[s].start:subproblems[s].finish) for s=1:S, k in keys(subproblem_solutions[s]))
                sum( λ[s, k]*value(subproblem_solutions[s][k]["vars"][:generator_output][gen, t]) for s=1:S, k in keys(subproblem_solutions[s]), t=subproblems[s].start:subproblems[s].finish)
                <= finish*generators[gen]["max_cf"])
    # Startup cost constraints
    @constraint(rmp, master_startup[gen in keys(generators), s=1:S-1], - sum(value.(subproblem_solutions[s][k]["vars"][:generator_on][gen, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    + sum(value.(subproblem_solutions[s+1][k]["vars"][:generator_on][gen, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1]))
    <=
    λ_startup[(s, s+1), gen])

    # Storage level constraints
    @constraint(rmp, master_storage[stg in keys(storage), s=1:S-1],
    sum(storage[stg]["efficiency"]*storage[stg]["capacity"]*value.(subproblem_solutions[s][k]["vars"][:storage_in][stg, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    - sum(storage[stg]["capacity"]*value.(subproblem_solutions[s][k]["vars"][:storage_out][stg, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    + sum(storage[stg]["energy"]*value.(subproblem_solutions[s][k]["vars"][:storage_level][stg, subproblems[s].finish])*λ[s, k] for k in keys(subproblem_solutions[s]))
    ==
    sum(storage[stg]["energy"]*value.(subproblem_solutions[s+1][k]["vars"][:storage_level][stg, subproblems[s+1].start])*λ[s+1, k] for k in keys(subproblem_solutions[s+1]))
    )

    master.model = rmp
    return true
end

function get_solution(master::Master_problem)
    subproblem_solutions = master.subproblem_solutions
    λ_values = value.(master.model.obj_dict[:λ]).data
    variables = Dict()
    for var in variable_names
        var_dict = Dict()
        for (s, sub) in subproblem_solutions
            var_dict[s] = 0*value.(subproblem_solutions[s][1]["vars"][var])
            for k in keys(sub)
                var_dict[s].data .+= λ_values[(s,k)]*value.(subproblem_solutions[s][k]["vars"][var]).data
            end
        end
        variables[var] = var_dict
    end
    return variables
end

function dual(master::Master_problem)
    # constraints = [name for (name, value) in master.model.obj_dict if typeof(value) <: JuMP.Containers.DenseAxisArray]
    constraints = [:max_cf, :master_startup, :master_storage]
    mu = Dict()
    for con in constraints
        mu[con] = Dict()
        for gen in master.model.obj_dict[con].axes[1]
            if con == :max_cf
                mu[con][gen] = dual.(master.model.obj_dict[con])[gen]
            elseif con in [:master_startup, :master_storage]
                mu[con][gen] = dual.(master.model.obj_dict[con])[gen,:].data
            end
        end
    end
    return mu
end

function create_solution!(master::Master_problem)
    master_solutions = get_solution(master)
    num_subproblems = master.S
    master.solution = Dict()
    for var in variable_names
        new_data = hcat([master_solutions[var][i].data for i in 1:num_subproblems]...) .|>
                        x->round(x; digits = 10) |> x->maximum([x,0])

        gens = master_solutions[var][1].axes[1]
        intervals = master.start:master.finish
        array = JuMP.Containers.DenseAxisArray(new_data, gens, intervals)
        master.solution[var] = array
    end
    return true
end

function graph_master(master::Master_problem)
    function create_demand_df(master)
        d = DataFrames.DataFrame(region = String[], interval = Float64[], value = Float64[])
        for (r, values) in master.demand
            for t in master.start:master.finish
                push!(d, [r t values[t]])
            end
            push!(d, [r master.finish+1 values[master.finish]])
        end
        return d
    end

    function create_variables_dict(master)
        generator_output = master.solution[:generator_output]
        storage_in = master.solution[:storage_in]
        storage_out = master.solution[:storage_out]
        regions = master.regions
        generation_vars = Dict()
        for (gen_name, gen_info) in master.inputs.generators
            generation_vars[gen_name] = Dict()
            generation_vars[gen_name]["generation"] = Dict()
            for t=master.start:master.finish
                generation_vars[gen_name]["generation"][t] = generator_output[gen_name,t]
            end
        end
        df = DataFrames.DataFrame(region = String[], name = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        for (gen_name, gen_info) in master.inputs.generators
            r = gen_info["region"]
            for t=master.start:master.finish
                push!(df, [r gen_name gen_info["capacity"] t max(gen_info["capacity"]*generation_vars[gen_name]["generation"][t],0)])
            end
        end
        storage_vars = Dict()
        for (stg_name, stg_info) in master.inputs.storage
            storage_vars[stg_name] = Dict()
            storage_vars[stg_name]["storage_in"] = Dict()
            storage_vars[stg_name]["storage_out"] = Dict()
            for t=master.start:master.finish
                storage_vars[stg_name]["storage_in"][t] = storage_in[stg_name,t]
                storage_vars[stg_name]["storage_out"][t] = storage_out[stg_name,t]
            end
        end
        stg_df = DataFrames.DataFrame(region = String[], name = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        # stg_df = DataFrames.DataFrame(region = String[], name = String[], out_or_in = String[], max_capacity = Float64[], interval = Int64[], generation = Float64[])
        for (stg_name, stg_info) in master.inputs.storage
            r = stg_info["region"]
            for t=master.start:master.finish
                # push!(stg_df, [r stg_name "in" stg_info["capacity"] t -stg_info["capacity"]*storage_vars[stg_name]["storage_in"][t]])
                # push!(stg_df, [r stg_name "out" stg_info["capacity"] t stg_info["capacity"]*storage_vars[stg_name]["storage_out"][t]])
                push!(stg_df, [r string(stg_name, "_in") stg_info["capacity"] t max(stg_info["capacity"]*storage_vars[stg_name]["storage_in"][t],0)])
                push!(stg_df, [r string(stg_name, "_out") stg_info["capacity"] t max(stg_info["capacity"]*storage_vars[stg_name]["storage_out"][t], 0)])
            end
        end
        stg_df = combine(groupby(stg_df, [:region, :name, :max_capacity, :interval]), [:generation => sum])
        stg_df = rename(stg_df, [:generation_sum => :generation])
    return [df; stg_df]
    end

    df = create_variables_dict(master)
    demand_df = create_demand_df(master)
    demand_df[:interval] = demand_df[:interval] .- 0.5

    plot(df, ygroup=:region,
       Geom.subplot_grid(
       layer(sort!(demand_df, rev=true), x=:interval,y=:value, ygroup=:region,Geom.step, Theme(default_color="black")),
       layer(sort!(df, rev=true),x=:interval,y=:generation,color=:name,ygroup=:region, Geom.bar)
       ))
end

function integer_heuristic(master::Master_problem)
    solution = master.solution
    return integer_heuristic(solution)
end

function integer_heuristic(node::Node)
    solution = node.solution
    return integer_heuristic(solution)
end

function integer_heuristic(solution)
    # First make generator_on variables binary
    axes = solution[:generator_output].axes
    integer_solution = Dict{Symbol, JuMP.Containers.DenseAxisArray{Float64, 2}}(:generator_output => solution[:generator_output])
    integer_generator_on = ceil.(solution[:generator_on])
    integer_solution[:generator_on] = integer_generator_on
    num_gens = size(integer_generator_on)[1]

    # Then calculate the new generator_startup values
    integer_generator_startup = ([integer_generator_on zeros(num_gens,1)] - [zeros(num_gens,1) integer_generator_on])[:, 1:end-1]
    integer_solution[:generator_startup] = max.(JuMP.Containers.DenseAxisArray(integer_generator_startup, axes...),0)

    integer_solution[:deficit] = solution[:deficit]
    integer_solution[:surplus] = solution[:surplus]

    generators = inputs["generators"]
    if sum(!(gen_info["capacity"]*integer_solution[:generator_output][gen, t] >= gen_info["mingen"]*integer_solution[:generator_on][gen, t]) for t in 1:T, (gen, gen_info) in generators) != 0
        # println("Integer solution does not satisfy the minimum generation constraint")
        return false
    end
    return integer_solution
end

function find_integer_objective_value(master::Master_problem)
    integer_solution = integer_heuristic(master)

    if integer_solution == false
        return false
    end

    start = master.start
    finish = master.finish
    generators = master.inputs.generators
    obj_value = 0
    obj_value += sum(gen_info["cost"]*gen_info["capacity"]*integer_solution[:generator_output][gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += sum(gen_info["startup"]*integer_solution[:generator_startup][gen, t] for (gen, gen_info) in generators, t in start:finish if t > 1)
    obj_value += 14700*sum(integer_solution[:deficit])
    obj_value += 1000*sum(integer_solution[:surplus])
    return obj_value
end

function find_integer_objective_value(node::Node)
    integer_solution = integer_heuristic(node)

    if integer_solution == false
        return false
    end

    start = 1
    finish = T
    generators = inputs["generators"]
    obj_value = 0
    obj_value += sum(gen_info["cost"]*gen_info["capacity"]*integer_solution[:generator_output][gen, t] for (gen, gen_info) in generators, t in 1:T)
    obj_value += sum(gen_info["startup"]*integer_solution[:generator_startup][gen, t] for (gen, gen_info) in generators, t in 2:T)
    obj_value += 14700*sum(integer_solution[:deficit])
    obj_value += 1000*sum(integer_solution[:surplus])
    return obj_value
end

function check_integer!(node::Node)
    if isnothing(node.solution)
        println("Node has not been solved yet")
        return false
    else
        n_iinf = sum(1 .- ((node.solution[:generator_on].data .== 0) .| (node.solution[:generator_on].data .== 1)))
        node.n_iinf = n_iinf
        if n_iinf == 0
            node.is_integer = true
            return true
        else
            return false
        end
    end
end

function most_fractional_index(master::Master_problem)
    length = UInt16(master.finish/master.S)
    ind = argmin(abs.(master.solution[:generator_on].data.-0.5))
    subproblem_id = div(ind.I[2], length)
    ind = CartesianIndex(ind.I[1], mod(ind.I[2], length))
    return (subproblem_id, ind)
end
function most_fractional_index(solution)
# function most_fractional_index(solution <: Dict)
    length = UInt16(T/S)
    ind = argmin(abs.(solution[:generator_on].data.-0.5))
    # Reduce the interval by 1 and then add 1 to get correct subproblem
    subproblem_id = div(ind.I[2]-1, length)+1
    ind = CartesianIndex(ind.I[1], mod(ind.I[2]-1, length)+1)
    return (subproblem_id, ind)
end


function branch(master::Master_problem, node::Node, branching_variables::Dict, node_counter)
    if node.is_integer
        return []
    end
    subproblem_id, ind = most_fractional_index(node.solution)
    variable_key = (subproblem_id, ind)
    if !((subproblem_id, ind) in keys(branching_variables))
        branching_variables[variable_key[1]][variable_key[2], 0] = findall(x -> value(x["vars"][:generator_on][ind]) != 0,
                                                       master.subproblem_solutions[subproblem_id])
        branching_variables[variable_key[1]][variable_key[2], 1] = findall(x -> value(x["vars"][:generator_on][ind]) != 1,
                                                       master.subproblem_solutions[subproblem_id])
    end
    left_node = Node(master, node_counter+1, node, variable_key, 0)
    right_node = Node(master, node_counter+2, node, variable_key, 1)
    return left_node, right_node
end

function constrain_master!(master::Master_problem, node::Node, branching_variables::Dict)
    for (s, fix_values) in enumerate(node.variable_values)
        for (variable_key, fix_value) in fix_values
            lambda_inds = branching_variables[variable_key[1]][variable_key[2], fix_value]
            # lambda_inds = branching_variables[variable_key, fix_value]
            for ind in lambda_inds
                fix(master.model[:λ][s, ind], 0; force = true)
            end
        end
    end
end

function master_heuristic(master::Master_problem)
    for x in master.model[:λ]
        set_binary(x)
    end
    optimize!(master.model)
    println("\n", objective_value(master.model))
    for x in master.model[:λ]
        unset_binary(x)
    end
    # println(value.(master.model[:λ]))
end


function sort_most_integer!(nodes::Array)
    sort!(nodes, by = x -> x.n_iinf, rev = true)
end

function sort_depth_first!(nodes::Array)
    sort!(nodes, by = x -> x.depth)
end

function sort_breadth_first!(nodes::Array)
    sort!(nodes, by = x -> x.depth, rev = true)
end

function print_node(tree_meta::Tree_meta)
    node_str = lpad(string(tree_meta.iteration_counter), 4)
    n_iinf_str = lpad(string(tree_meta.n_iinf), 4)
    nodes_left_str = lpad(string(tree_meta.n_nodes_left), 4)
    lower_bound_str = lpad(string(round(tree_meta.lower_bound, digits = 2)), 2)
    mip_gap = round(tree_meta.mip_gap, digits = 4)*100
    print(node_str, "     ", n_iinf_str, "           ", nodes_left_str, "        ", mip_gap, "     ", lower_bound_str)
    if tree_meta.integer_found
        print("   ", lpad(string(round(tree_meta.best_integer, digits = 4)), 1))
    end
end

function find_node_relaxation(node::Node)
    start = 1
    finish = T
    generators = inputs["generators"]
    solution = node.solution
    obj_value = 0
    obj_value += sum(gen_info["cost"]*gen_info["capacity"]*solution[:generator_output][gen, t] for (gen, gen_info) in generators, t in start:finish)
    obj_value += sum(gen_info["startup"]*solution[:generator_startup][gen, t] for (gen, gen_info) in generators, t in start:finish if t > 1)
    obj_value += 14700*sum(solution[:deficit])
    obj_value += 1000*sum(solution[:surplus])
    return obj_value
end

function deep_dive(node::Node, subproblems)
    length = Int64(T/S)
    fixed_variables = round.(node.solution[:generator_on]).data
    for (s, subproblem) in enumerate(subproblems)
        for ind in CartesianIndices(fixed_variables)
            modified_ind = CartesianIndex(ind.I[1], mod(ind.I[2]-1, length)+1)
            unset_binary(subproblem.model[:generator_on][modified_ind])
            fix(subproblem.model[:generator_on][modified_ind], fixed_variables[ind])
        end
    end
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
    for (s, subproblem) in enumerate(subproblems)
        for ind in CartesianIndices(fixed_variables)
            modified_ind = CartesianIndex(ind.I[1], mod(ind.I[2]-1, length)+1)
            unfix(subproblem.model[:generator_on][modified_ind], fixed_variables[ind])
            set_binary(subproblem.model[:generator_on][modified_ind])
        end
    end

end

end
