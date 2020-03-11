
# Solving the unit commitment problem
# using CPLEX
using GLPK
using JuMP
using YAML
using DataFrames
using Gadfly
using Random

GLPK.jl_set_preemptive_check(false)

# First solve D_i
Random.seed!(100)

T = 3*48
# T = 50*48
generators = YAML.load(open("generators-regions.yml"))
regions = keys(generators) |> collect
demand = Dict()
directions = ["import", "export"]
for t=1:T
    rand_addition = rand()
    for r in regions
        if t==1
            demand[r] = Float64[]
        end
        append!(demand[r], 45+40*(rand()-0.5)*rand_addition)
    end
end

demand["NSW"][Int(T/2)+1] = demand["NSW"][Int(T/2)]-30

# original = JuMP.Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_PREIND = 0, CPX_PARAM_PREDUAL = -1, CPX_PARAM_MIPDISPLAY = 2))
original = JuMP.Model(with_optimizer(GLPK.Optimizer, msg_lev = GLPK.MSG_ALL))

# Generator variables
@variable(original, 1>=generator_output[r=regions, keys(generators[r]["generators"]), t=1:T]>=0)
@variable(original, generator_on[r=regions, keys(generators[r]["generators"]), t=1:T], Bin)
@variable(original, generator_startup[r=regions, keys(generators[r]["generators"]), t=1:T], Bin)

# Interconnector variables
@variable(original, 1>=interconnector[r=regions, keys(generators[r]["interconnectors"]), t=1:T, direction = directions]>=0)
# @variable(original, 1>=interconnector_import[r=regions, keys(generators[r]["interconnectors"]), t=1:T]>=0)

@objective(original, Min, sum(d["cost"] * d["capacity"] * generator_output[r, gen, t] + d["startup"]*generator_startup[r, gen,t] for r in regions, (gen,d) in generators[r]["generators"],t=1:T))

# Demand constraints
@constraint(original, demand_constraint[r=regions,t=1:T], sum(d["capacity"]*generator_output[r,gen,t] for (gen,d) in generators[r]["generators"]) -
    sum(i["capacity"]*interconnector[r, int, t, "export"] for (int, i) in generators[r]["interconnectors"]) + sum(i["capacity"]*interconnector[r, int, t, "import"] for (int, i) in generators[r]["interconnectors"]) == demand[r][t])
# @constraint(original, demand_constraint[r=regions,t=1:T], sum(d["capacity"]*generator_output[r,gen,t] for (gen,d) in generators[r]["generators"]) == demand[r][t])

# Min gen constraints
@constraint(original, min_gen[r=regions, gen in keys(generators[r]["generators"]), t=1:T], generators[r]["generators"][gen]["capacity"]*generator_output[r, gen, t] >= generators[r]["generators"][gen]["mingen"]*generator_on[r, gen, t])
# On variable constraint
@constraint(original, on_contraint[r=regions, gen in keys(generators[r]["generators"]), t=1:T], generator_output[r, gen, t] <= generator_on[r, gen, t])
# Ramp rate constraints
@constraint(original, ramp_down[r=regions, gen in keys(generators[r]["generators"]), t=1:T-1], generators[r]["generators"][gen]["capacity"]*(generator_output[r, gen, t] - generator_output[r, gen, t+1]) <= generators[r]["generators"][gen]["ramp"])
@constraint(original, ramp_up[r=regions, gen in keys(generators[r]["generators"]), t=1:T-1], generators[r]["generators"][gen]["capacity"]*(generator_output[r, gen, t] - generator_output[r, gen, t+1]) >= -generators[r]["generators"][gen]["ramp"])
# Start up constraint
@constraint(original, start_up[r=regions, gen in keys(generators[r]["generators"]), t=2:T], generator_startup[r, gen,t]>=generator_on[r, gen,t]-generator_on[r, gen,t-1])

# Interconnector loss constraints
@constraint(original, interconnector_loss[r=regions, int in keys(generators[r]["interconnectors"]), t=1:T], (1-generators[r]["interconnectors"][int]["loss"])*interconnector[r, int, t, "export"] == interconnector[generators[r]["interconnectors"][int]["to"], int, t, "import"])

optimize!(original)
# solve(original)

function create_variables_dict(generators)
    generation_vars = Dict()
    for r in regions
        generation_vars[r] = Dict()
        for gen in keys(generators[r]["generators"])
            generation_vars[r][gen] = Dict()
            generation_vars[r][gen]["generation"] = Dict()
            generation_vars[r][gen]["on"] = Dict()
            for t=1:T
                generation_vars[r][gen]["generation"][t] = value.(generator_output[r, gen,t])
                generation_vars[r][gen]["on"][t] = value.(generator_on[r, gen,t])
            end
        end
    end
    return generation_vars
end

generation_vars = create_variables_dict(generators)

# Create the dataframe
function create_generation_df(generators)
    df = DataFrames.DataFrame(region = String[], generator_name = String[], cost = Float64[], max_capacity = Float64[], interval = Int64[], generation = Float64[], demand = Float64[])
    for r in regions
        for (gen, d) in generators[r]["generators"]
            for t=1:T
                push!(df, [r gen d["cost"] d["capacity"] t d["capacity"]*generation_vars[r][gen]["generation"][t] demand[r][t]])
            end
        end
    end
    return df
end

df = create_generation_df(generators)

df[:generation] = max.(df[:generation],0)

function create_demand_df(demand)
    df = DataFrames.DataFrame(region = String[], interval = Float64[], value = Float64[])
    for (r, values) in demand
        for t in 1:T
            push!(df, [r t values[t]])
        end
        push!(df, [r T+1 values[T]])
    end
    return df
end

function get_interconnectors(interconnector)
    df = DataFrames.DataFrame(region = String[], name = String[], interval = Float64[], direction = String[], value = Float64[])
    for r in regions
        for int in keys(generators[r]["interconnectors"])
            for t in 1:T
                for direction in directions
                    v = value.(interconnector)[r, int, t, direction]
                    push!(df, [r int t direction v])
                end
            end
        end
    end
    return df
end

demand_df = create_demand_df(demand)
demand_df[:interval] = demand_df[:interval] .- 0.5
plot(df, ygroup=:region,
   Geom.subplot_grid(
   layer(sort!(demand_df, rev=true), x=:interval,y=:value, ygroup=:region,Geom.step, Theme(default_color="black")),
   layer(sort!(df, rev=true),x=:interval,y=:generation,color=:generator_name,ygroup=:region, Geom.bar))
   )
# plot(df,
#    layer(demand_df, x=:interval,y=:value, ygroup=:region,Geom.subplot_grid(Geom.step), Theme(default_color="black")),
#    layer(sort!(df, rev=true),x=:interval,y=:generation,color=:generator_name,ygroup=:region, Geom.subplot_grid(Geom.bar)))
# plot(df,
#    layer(demand_df, generator_output=:interval,generator_on=:demand,Geom.step, Theme(default_color="black")),
#    layer(sort!(df, rev=true),generator_output=:interval,generator_on=:generation,color=:generator_output,Geom.bar))
