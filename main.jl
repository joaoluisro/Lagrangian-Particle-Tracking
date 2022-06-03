using NetCDF

function ReadField(filename :: String)

    u_vel = ncread(filename, "uo")
    #v_vel = ncread("dataset-baltic-daily.nc", "vo")
    return u_vel
end


include("particle.jl")
v1(x,y) = [x, y]/(x*x + y*y)
v2(x,y,t) = [-y,x]/(y*y + x*x)
particle_initial = rand(Uniform(-15,15),500, 2)
domain_range = (-15:60,-50:50)
n_iter = 300

result = Advect(particle_initial,
                v2,
                domain_range,
                time_step = 0.01,
                n_iter = 300)

Animate(n_iter,
        domain_range,
        result)
