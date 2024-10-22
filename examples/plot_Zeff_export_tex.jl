
using ADAS, Plots, OrderedCollections, Format, PGFPlotsX

impurities = [:Ne, :Ar, :Cr, :Kr, :W]
color_ = [:blue, :green, :red, :black, :cyan, :red]
Lz = OrderedDict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
Zeff = Dict(imp => ADAS.get_Zeff(imp) for imp in impurities)
pgfplotsx()
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle=:box, palette=:darktest, size=(400, 400))
for (i, Lz_) in enumerate(values(Lz))
    plot!(Lz_, Te=Log10Range(0.5, 4.5, 1000), xlim=[5.0, 2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i])
end
# y = Log10Range(-34, -31, 4)
# x = Log10Range(1, 4, 4)
# yticks!(p, y)
# xticks!(p, x)
# plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)
# #%%
# using NumericalIntegration
# Prad_tot = Dict(k => NumericalIntegration.integrate(r_omp, v) for (k, v) in Pvol_rad)
# fraction_imp = 0.0:0.001:0.005
# Zeff_ped = Dict(imp => [ADAS.get_Zeff(imp)(f, ne_ped, Te_ped) for f in fraction_imp] for imp in impurities)
using LaTeXStrings
p = plot()
Te = Log10Range(1, 4, 100)
ne = 1e19
fraction = [0.001, 0.01]
color_ = [:blue, :green, :red, :black, :cyan, :red]
linestyle_ = [:solid,:dash]
plot!(xlabel="Te[eV]", ylabel=L"$Z_{eff} = \frac{\sum n_a Z_a^2}{n_e}")
plot!(title="Zeff in C-R equilibrium\n", titlefontsize=8)
for (i,imp) in enumerate([:Ne,:Ar,:Kr])
    for (j,f) in enumerate(fraction)
        plot!(Te, [Zeff[imp](f, ne, Te_) for Te_ in Te], label="$imp : $(f*100)%", lw=1.0, color=color_[i],linestyle=linestyle_[j])
    end
end
fs =12
plot!(framestyle=:box, ylims=(0, 5.0), xtickfontsize=fs, ytickfontsize=fs, xguidefontsize=fs, yguidefontsize=fs, legendfontsize=fs)
#save("zeff.tikz", p)