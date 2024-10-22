
using ADAS, Plots, OrderedCollections, Format
gr()

using PyCall
aurora = pyimport("aurora")
np =pyimport("numpy")
Te_eV = np.logspace(np.log10(1), np.log10(1e4), 1000)
ne_cm3 = 10e13 * np.ones_like(Te_eV)
imp = "Kr"
line_rad_tot, cont_rad_tot = aurora.get_cooling_factors(imp, ne_cm3, Te_eV, plot=true, ion_resolved=false)
files = ["scd", "acd"]
atom_data_eq = aurora.get_atom_data(imp, files)
_Te, fz = aurora.get_frac_abundances(atom_data_eq, ne_cm3, Te_eV)
#=
import aurora
import numpy as np
import matplotlib.pyplot as plt

plt.ion()


# scan Te and fix a value of ne
Te_eV = np.logspace(np.log10(1), np.log10(1e4), 1000)
ne_cm3 = 10e13 * np.ones_like(Te_eV)

imp = "Kr"

# basic cooling curve, considering ionization equilibrium between charge states
line_rad_tot, cont_rad_tot = aurora.get_cooling_factors(
    imp, ne_cm3, Te_eV, plot=True, ion_resolved=False
)


aurora.get_frac_abundances(imp, ne_cm3, Te_eV, plot=True)
files = ["scd", "acd"]
atom_data_eq = aurora.get_atom_data(imp, files)


_Te, fz = aurora.get_frac_abundances(atom_data_eq, ne_cm3, Te_eV)
=#
impurities = [:Kr]
color_ = [:orange, :gray, :blue, :black, :cyan, :red]

Lz = OrderedDict(imp => ADAS.get_cooling_rates(imp; plt_year="41") for imp in impurities);
using Plots
Log10Range(args...; kw...) = 10 .^ LinRange(args...; kw...)
p = plot(layout=1, framestyle=:box, palette=:darktest, size=(400, 400))
for (i, Lz_) in enumerate(values(Lz))
    plot!(Lz_, Te=Log10Range(0, 4, 1000), xlim=[5.0, 2e4], ylim=[1e-34, 5e-31], xscale=:log10, yscale=:log10, layout=1, color=color_[i])
end
plot!(Te_aurora, Lz_aurora, label="aurora Kr")
plot!()
y = Log10Range(-34, -31, 4)
x = Log10Range(1, 4, 4)
yticks!(p, y)
xticks!(p, x)
plot!(legendposition=:bottomright, legendfontsize=12, legendcolumns=2)