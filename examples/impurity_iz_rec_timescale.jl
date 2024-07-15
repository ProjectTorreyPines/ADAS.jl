#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 impurity_iz_rec_timescale.jl (c) 2024=#

using ADAS, Plots
pythonplot()
impurities = [:C, :Ne, :Si, :Ar, :Kr, :Xe, :W]
Lz = [ADAS.get_cooling_rates(imp) for imp in impurities];
Rrad = [ADAS.get_radiation_rates(imp) for imp in impurities];
Zeff = [ADAS.get_Zeff(imp) for imp in impurities];
Zmean = Dict(imp => ADAS.get_Zmean(imp) for imp in impurities);
Rrec = Dict(imp => ADAS.get_recombination_rate(imp) for imp in impurities);
Riz = Dict(imp => ADAS.get_ionization_rate(imp) for imp in impurities)
# Zmean = Dict(Lz_.imp => [ADAS.get_Zmean(Lz_.imp)(ne_, Te_) for (Te_, ne_) in zip(Te, ne)] for Lz_ in Lz)
using Format
using LaTeXStrings
Te__ = 10 .^ LinRange(0.0, 5.0, 100)
ne__ = 10 .^ LinRange(19, 22, 100)
Te, ne = meshgrid(Te__, ne__)
#plot(layout=(length(Zmean),2))
τ_iz = Dict(imp => [1.0 / (Riz_.rate((Zmean[imp](ne_, Te_)), ne_, Te_) * ne_) for (Te_, ne_) in zip(Te, ne)] for ((imp, Riz_)) in Riz)
τ_rec = Dict(imp => [1.0 / (Rrec_.rate((Zmean[imp](ne_, Te_)), ne_, Te_) * ne_) for (Te_, ne_) in zip(Te, ne)] for ((imp, Rrec_)) in Rrec)
τ_max = Dict(imp => map(x -> maximum(x), zip(τ_iz[imp], τ_rec[imp])) for imp in impurities)
# plot(layout=length(impurities), size=(1200, 800))
# contourf!(ne__, Te__, log10.(τ_rec[:Kr]), color=:turbo, clabels=true, cbar=true, subplot=1, xscale=:log10, yscale=:log10, clims=(-6, 0), aspect_ratio=:auto)
# contourf!(ne__, Te__, log10.(τ_iz[:Kr]), color=:turbo, clabels=true, cbar=true, subplot=2, xscale=:log10, yscale=:log10, clims=(-6, 0))
# contourf!(ne__, Te__, (τ_iz[:Kr] ./ τ_rec[:Kr]), color=:turbo, clabels=true, cbar=true, subplot=3, xscale=:log10, yscale=:log10, clims=(0.9, 1.1))
# plot!()
plot(; layout=length(impurities), size=(1200, 800))
colorbar_ticks = (collect(-6:1:1.0), format.("{:2.2e}", 10 .^ collect(-6:1:1.0)))
for (i, imp) in enumerate(impurities)
    contourf!(
        ne__,
        Te__,
        log10.(τ_max[imp]);
        annot=1e-1,
        levels=collect(-6:1),
        color=:turbo,
        xlabel=L"$n_e [m^{-3}]$",
        ylabel=L"$T_e [eV]$",
        clabels=false,
        cbar=true,
        subplot=i,
        xscale=:log10,
        colorbar_ticks=colorbar_ticks,
        yscale=:log10,
        title=L"$\max(\tau_{iz},\tau_{rec})[s]$" * ": $imp",
        clim=(-6, 1.0)
    )
    plot!(; framestyle=:box, grid=true, gridlinewidth=3, subplot=i)

    vline!(10 .^ LinRange(19, 22, 4); color="gray", label="", style=:dash, subplot=i)
    hline!(10 .^ LinRange(0, 5, 6); color="gray", label="", style=:dash, subplot=i)
    fontsize = 12
    annotate!(1e19, 50000, text("$imp", fontsize, :Courier, :left); subplot=i)
end
plot!()

# plot(title="Average charge s