# Solving the unit commitment problem
# using CPLEX
using GLPK
using JuMP
using YAML
using DataFrames
using Gadfly
using Random


# include("Decomposition.jl")
# using .Decomposition
Random.seed!(100)

function create_master_problem(subproblem_solution, subproblem_objective_values, voll_cost = voll_cost)
    S = length(subproblem_solution)
    # master = JuMP.Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_PREIND = 0, CPX_PARAM_PREDUAL = -1, CPX_PARAM_MIPDISPLAY = 1))
    master = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))

    @variable(master, 1>=lambda[k=1:K, s=1:S]>=0)
    @variable(master, 1>=lambda_voll[k=1:K]>=0)

    # Objective function
    @objective(master, Min, sum(subproblem_objective_values[s][k] * lambda[k,s] for k=1:K, s=1:S) +
                            sum(voll_cost[k]*lambda_voll[k] for k=1:K))

    @constraint(master, ramp_down[r=regions, gen in keys(generators[r]["generators"]), k=1:K-1],
        generators[r]["generators"][gen]["capacity"] * sum(subproblem_solution[s][k+1][r][gen][1]*lambda[k+1,s] for s=1:S) >=
        generators[r]["generators"][gen]["capacity"] * sum(subproblem_solution[s][k][r][gen][T]*lambda[k,s] for s=1:S) - generators[r]["generators"][gen]["ramp"])
    @constraint(master, ramp_up[r=regions, gen in keys(generators[r]["generators"]), k=1:K-1],
        generators[r]["generators"][gen]["capacity"] * sum(subproblem_solution[s][k+1][r][gen][1]*lambda[k+1,s] for s=1:S) <=
        generators[r]["generators"][gen]["capacity"] * sum(subproblem_solution[s][k][r][gen][T]*lambda[k,s] for s=1:S) + generators[r]["generators"][gen]["ramp"])

    @constraint(master, convexity_constraint[k=1:K], sum(lambda[k,s] for s=1:S) + sum(lambda_voll[k])== 1)
    optimize!(master)

    duals = Dict()
    # duals["ramp_up"] = getdual.(ramp_up)
    # duals["ramp_down"] = getdual.(ramp_down)
    duals["ramp_up"] = JuMP.dual.(ramp_up)
    duals["ramp_down"] = JuMP.dual.(ramp_down)
    lambda = getvalue.(lambda)
    convexity_duals = JuMP.dual.(convexity_constraint)
    println(convexity_duals)
    return master, duals, lambda, convexity_duals
end

T = 48
K = 3
generators = YAML.load(open("generators-regions.yml"))
regions = keys(generators) |> collect
demand = Dict()
directions = ["import", "export"]
for t=1:(K*T)
    rand_addition = rand()
    for r in regions
        if t==1
            demand[r] = Float64[]
        end
        append!(demand[r], 45+40*(rand()-0.5)*rand_addition)
    end
end

for r in keys(demand)
    demand[r] = reshape(demand[r], T, K)
end

demand["NSW"][T+1] = demand["NSW"][T]-30

# Set up the master problem
joining_sections = []
for k=1:K-1
    push!(joining_sections, (k,k+1))
end

mu = Dict()
for k=1:K-1
    mu[(k,k+1)] = Dict()
    for r in regions
        mu[(k,k+1)][r] = Dict()
        for gen in keys(generators[r]["generators"])
            mu[(k,k+1)][r][gen] = Dict("ramp_up" => 0, "ramp_down" => 0)
            # mu[(k,k+1)][r][gen]["ramp_up"] = 0
            # mu[(k,k+1)][r][gen]["ramp_down"] = 0
        end
    end
end
# mu = Dict(
# "ramp_up" => Dict(),
# "ramp_down" => Dict(),
# )
# for (key, value) in mu
#     for k=1:K-1
#         value[(k,k+1)] = Dict()
#         for r in regions
#             value[(k,k+1)][r] = Dict()
#             for gen in keys(generators[r]["generators"])
#                 value[(k,k+1)][r][gen] = 0
#             end
#         end
#     end
# end

## Create the subproblems
subproblems = []
subproblem_solution = []
subproblem_objective_values = []
convexity_duals = zeros(K)
S = 1

voll_cost = Array{Float64}(undef, K)
for k=1:K
    voll_cost[k] = sum.(v[:,k]*generators[region]["generators"]["dsm_voll"]["cost"] for (region, v) in demand) |> sum
end
# master, mu, lambda, convexity_duals = create_master_problem(subproblem_solution, subproblem_objective_values)

function create_variables_df(subproblem_solutions)
    df = DataFrames.DataFrame(region = String[], generator_name = String[], interval = Int64[], generation = Float64[])
    for (k, v) in enumerate(subproblem_solutions)
        for (r, gens) in v
            for (gen, output) in gens
                for t=1:T
                    push!(df, [r gen ((k-1)*T + t) generators[r]["generators"][gen]["capacity"]*subproblem_solutions[k][r][gen][t]])
                end
            end
        end
    end
    return df
end

master, mu, lambda, convexity_duals = create_master_problem(subproblem_solution, subproblem_objective_values)

mu = Dict()
for k=1:K-1
    mu[(k,k+1)] = Dict()
    for r in regions
        mu[(k,k+1)][r] = Dict()
        for gen in keys(generators[r]["generators"])
            mu[(k,k+1)][r][gen] = Dict("ramp_up" => 0, "ramp_down" => 0)
            # mu[(k,k+1)][r][gen]["ramp_up"] = 0
            # mu[(k,k+1)][r][gen]["ramp_down"] = 0
        end
    end
end
# subproblem1 = Decomposition.subproblem(1, T, demand, generators)
# subproblem_solution, subproblem_objective_values = Decomposition.solve_subproblem(subproblem1, mu[(1,2)], convexity_dual[1])

# subproblem_solution, subproblem_objective_values = Decomposition.solve_subproblem(subproblem_solution, mu)
#
# master, mu = create_master_problem(subproblem_solution, subproblem_objective_values)
