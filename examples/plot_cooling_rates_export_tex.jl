
using ADAS, Plots, OrderedCollections, Format, PGFPlotsX
gr()
impurities = [:Ne, :Ar, :Cr, :Kr, :W]
color_ = [:orange, :gray, :blue, :black, :cyan, :red]
Lz = OrderedDict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
Zeff = Dict(imp => ADAS.get_Zeff(imp) for imp in impurities)
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle=:box, palette=:darktest, size=(400, 400))
for (i, Lz_) in enumerate(values(Lz))
    plot!(Lz_, Te=Log10Range(0.5, 4.5, 1000), xlim=[5.0, 2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i])
end
y = Log10Range(-34, -31, 4)
x = Log10Range(1, 4, 4)
yticks!(p, y)
xticks!(p, x)
plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)
#%%
using NumericalIntegration
Prad_tot = Dict(k => NumericalIntegration.integrate(r_omp, v) for (k, v) in Pvol_rad)
fraction_imp = 0.0:0.001:0.005
Zeff_ped = Dict(imp => [ADAS.get_Zeff(imp)(f, ne_ped, Te_ped) for f in fraction_imp] for imp in impurities)

Zeff_sep = Dict(imp => [ADAS.get_Zeff(imp)(f, ne_sep, Te_sep) for f in fraction_imp] for imp in impurities)
plot(layout=2)
plot!(title="Zeff at pedestal in C-R equilibrium\n (nₑ=$ne_ped,Tₑ=$Te_ped) ", subplot=1, titlefontsize=8)

for (k, v) in Zeff_ped
    plot!(Te, Zeff, label="$k", subplot=1)
end

plot!(xlabel="fraction of impurity", ylabel="Zeff = ∑nₐZₐ^2/nₑ", subplot=1, ylim=[0, 10])

plot!(title="Zeff at sep in C-R equilibrium\n (nₑ=$ne_sep,Tₑ=$Te_sep) ", subplot=2, titlefontsize=8)
Zeff_sep = Dict(imp => [ADAS.get_Zeff(imp)(f, ne_sep, Te_sep) for f in fraction_imp] for imp in impurities)
Zeff = Dict(imp => ADAS.get_Zeff(imp) for imp in impurities)
for (k, v) in Zeff_sep

    plot!(fraction_imp, v, label="$k", subplot=2)
end

plot!(xlabel="fraction of impurity", ylabel="Zeff = ∑nₐZₐ²/nₑ", subplot=2, ylim=[0, 10])