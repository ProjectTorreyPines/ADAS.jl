
using ADAS, Plots
gr()
impurities = [:Kr, :W, :Ar, :Ne, :C, :Si, :Cr]
Lz = Dict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
plot(layout=1, framestyle=:box,palette)
for Lz_ in values(Lz)
    plot!(Lz_, Te=Log10Range(0, 4, 1000), ylim=[1e-35, 5e-31], yscale=:log10, layout=1)
end
plot!()