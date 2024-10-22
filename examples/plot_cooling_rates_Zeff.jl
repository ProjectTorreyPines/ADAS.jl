
using ADAS, Plots, OrderedCollections, Format
plotly()
impurities = [:Ne, :Ar, :Cr, :Kr,:W]
color_ = [:blue, :green, :red, :black, :cyan]
Lz = OrderedDict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
Zeff = OrderedDict(imp => ADAS.get_Zeff(imp) for imp in impurities);
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle=:box, palette=:darktest, size=(400, 400))
for (i,Lz_) in enumerate(values(Lz))
    plot!(Lz_, Te=Log10Range(0.5, 4.5, 1000), xlim=[5.0,2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i])
end
Te = Log10Range(0.5, 4.5, 1000)

for (i, (Lz_,Zeff_)) in enumerate(zip(values(Lz),values(Zeff)))
    plot!([Zeff_(0.001, 1e20, Te_) for Te_ in Te], [Lz_.Lztot(1e20, Te_) for Te_ in Te], xlim=[1.0, 5], ylim=[1e-34, 5e-31], yscale=:log10, layout=1, color=color_[i], label=String(Lz_.imp))
end
y = Log10Range(-34, -31, 4)
x = Log10Range(1, 4, 4)
yticks!(p, y)
#xticks!(p, x)
plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)
#savefig("cooling_rates.tikz")