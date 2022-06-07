using NetCDF, Logging, Distributions, Plots, PlotlyJS

function Read2DFieldDeprecated(filename :: String;
                     field_names :: Tuple{String},
                     u_field :: String,
                     v_field :: String)

    if split(filename, ".")[end] != "nc"
        @error "Not a NETCDF file"
    end
    u_vel = ncread(filename, u_field)
    v_vel = ncread(filename, v_field)
    return (u_vel, v_vel)
end

function Read2DField(filename :: String;
                     variable_dict :: Dict{String, String})

    if split(filename, ".")[end] != "nc"
        @error "Not a NETCDF file"
    end

    for (variable, variable_name) in variable_dict
        try

        catch
            @error "Reading " variable
        end
    end
    u_vel = ncread(filename, u_field)
    v_vel = ncread(filename, v_field)
    return (u_vel, v_vel)
end

function Bilinear2DField(position :: Array{Float64},
                         fieldData :: Array{Float32, 2})
    x,y = position[1], position[2]
    grid_x, grid_y = ceil(Int64, x), ceil(Int64, y)

    v_grid = (fieldData[grid_x, grid_y],
              fieldData[grid_x + 1, grid_y],
              fieldData[grid_x, grid_y + 1],
              fieldData[grid_x + 1, grid_y + 1])

    return (v_grid[1] * (grid_x + 1 - x) * (grid_y + 1 - y)) +
           (v_grid[2] * (x - grid_x) * (grid_y + 1 - y)) +
           (v_grid[3] * (grid_x + 1 - x) * (y - grid_y)) +
           (v_grid[4] * (x - grid_x) * (y - grid_y))
end


function EulerUpdateBilinear(position :: Array{Float64},
                            u_field :: Array{Float32, 2},
                            v_field :: Array{Float32, 2},
                            time_step :: Float64,
                            domain_range :: Tuple{UnitRange{Int64}, UnitRange{Int64}},
                            threshold)

    x_range, y_range = collect(domain_range)
    # Checks boundary conditions
    in_range = true
    vx, vy = Bilinear2DField(position, u_field), Bilinear2DField(position, v_field)

    velocity = [vx, vy]
    valid_vel = minimum(velocity) > -10.0 && maximum(velocity) < 100.0
    if !valid_vel || !in_range
        @info "invalid" velocity[1], velocity[2]
        println(domain_range)
        println(position)
    end

    return position + (valid_vel && in_range)*(time_step * velocity)
end

function AdvectDataset(particle_initial :: Matrix{Float64},
                       u_vel :: Array{Float32, 3},
                       v_vel :: Array{Float32, 3},
                       domain_range :: Tuple{UnitRange{Int64}, UnitRange{Int64}};
                       time_step :: Float64 = 0.1,
                       missing_value = -Inf)


    N_PARTICLE = size(particle_initial)[1]
    # Initialize trajectory history

    n_iter = size(u_vel)[end]
    trajectory = zeros(Float64, N_PARTICLE, 2, n_iter)
    trajectory[:,:,1] = particle_initial

    for i in 2:n_iter
        updated = [EulerUpdateBilinear(collect(p),
                    u_vel[:,:,i-1],
                    v_vel[:,:,i-1],
                    time_step,
                    domain_range,
                    -500.00)
                    for p in eachrow(trajectory[:,:, i - 1])]
        trajectory[:,:,i] = mapreduce(permutedims, vcat, updated)
    end

    return trajectory
end

function Animate(n_iter :: Int64,
                domain_range :: Tuple{UnitRange{Int64}, UnitRange{Int64}},
                result :: Array{Float64, 3})

    x_range, y_range = collect(domain_range)


    animation = @animate for i in 1:n_iter
        state = result[:,:,i]
        Plots.scatter(xlims=(x_range[1], x_range[end]), ylims=(y_range[1], y_range[end]), state[:,1], state[:,2], legend=false)
    end

    return animation
end
