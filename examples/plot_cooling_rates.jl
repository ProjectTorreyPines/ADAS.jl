
using ADAS, Plots
impurities = [:Be, :Ne, :Ar, :Kr, :Cr, :W]
color_ = [:blue,:green, :red, :black,:cyan]
Lz = Dict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle=:box)
for Lz_ in values(Lz)
    plot!(Lz_, Te=Log10Range(0, 4, 1000), ylim=[1e-35, 5e-31], yscale=:log10, layout=1)
end
y = Log10Range(-34, -31, 4)
x = Log10Range(1, 4, 4)
yticks!(p, y)
xticks!(p, x)
plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)