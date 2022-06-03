using Distributions, Plots


function ForwardEulerUpdate(position :: Array{Float64},
                            velocity_field :: Function,
                            time_step :: Float64,
                            domain_range :: Tuple{UnitRange{Int64}, UnitRange{Int64}},
                            t :: Float64)

    x_range, y_range = collect(domain_range)
    # Checks boundary conditions
    in_range = (position[1] in x_range) && (position[2] in y_range)

    # Euler forward update:  p_n+1 = p_n + t * v(p, t)
    return position + in_range*(time_step * velocity_field(position[1], position[2], t))
end

function Advect(particle_initial :: Matrix{Float64},
                velocity_field :: Function,
                domain_range :: Tuple{UnitRange{Int64}, UnitRange{Int64}};
                time_step :: Float64 = 0.1,
                n_iter :: Integer = 100)

    N_PARTICLE = size(particle_initial)[1]
    # Initialize trajectory history
    trajectory = zeros(Float64, N_PARTICLE, 2, n_iter)
    trajectory[:,:,1] = particle_initial

    # Update position
    for i in 2:n_iter
        updated = [ForwardEulerUpdate(collect(p), velocity_field, time_step, domain_range, Float64(i))
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
        scatter(xlims=(x_range[1], x_range[end]), ylims=(y_range[1], y_range[end]), state[:,1], state[:,2], legend=false)
    end

    return animation
end
