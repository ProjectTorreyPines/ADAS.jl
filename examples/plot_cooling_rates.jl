
using ADAS, Plots, OrderedCollections, Format
pgfplotsx()
impurities = [:Ne, :Ar, :Kr]
color_ = [:blue,:green, :red, :black,:cyan]
Lz = OrderedDict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle = :box, palette=:darktest, size=(400, 400))
# for (i,Lz_) in enumerate(values(Lz))
#     plot!(Lz_, Te=Log10Range(0.5, 4.5, 1000), xlim=[5.0,2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i])
# end
Te=Log10Range(0.5, 4.5, 100)

for (i, Lz_) in enumerate(values(Lz))
    plot!(Te, [Lz_.Lztot(1e20, Te_) for Te_ in Te]./Te.^(7/4), xlim=[5.0, 2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i], label=String(Lz_.imp))
end
y = Log10Range(-34, -31, 4)
x = Log10Range(1, 4, 4)
yticks!(p, y)
xticks!(p, x)
plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)
savefig("cooling_rates.tikz")