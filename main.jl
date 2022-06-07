include("field.jl")
using PlotlyJS
x_initial = rand(Uniform(300, 400), 200)
y_initial = rand(Uniform(150, 200), 200)
particle_initial = hcat(x_initial, y_initial)

domain_range = (0:763, 0:774)
n_days = 31

filename = "dataset-baltic-daily.nc"
u, v = Read2DField(filename, "uo", "vo")
u = dropdims(u; dims = 3)
v = dropdims(v; dims = 3)

result = AdvectDataset(particle_initial,
                u,
                v,
                domain_range,
                time_step = 50.0)

anim = Animate(31, domain_range,result)
function viewTopography()

        layout = Layout(
                title="Topography",
                autosize=false,
                width=800,
                height=800,
                margin=attr(l=65, r=50, b=65, t=90))
return
end
