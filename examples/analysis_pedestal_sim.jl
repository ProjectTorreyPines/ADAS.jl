#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 impurity_iz_rec_timescale.jl (c) 2024=#

using ADAS, Plots
gr()
impurities = [:C, :Ne, :Si, :Ar, :Kr, :Xe, :W]
Lz = Dict(imp => ADAS.get_cooling_rates(imp) for imp in impurities);
Rrad = [ADAS.get_radiation_rates(imp) for imp in impurities];
Zeff = [ADAS.get_Zeff(imp) for imp in impurities];
Zmean = Dict(imp => ADAS.get_Zmean(imp) for imp in impurities);
Rrec = Dict(imp => ADAS.get_recombination_rate(imp) for imp in impurities);
Riz = Dict(imp => ADAS.get_ionization_rate(imp) for imp in impurities)
# Zmean = Dict(Lz_.imp => [ADAS.get_Zmean(Lz_.imp)(ne_, Te_) for (Te_, ne_) in zip(Te, ne)] for Lz_ in Lz)
using Format
using LaTeXStrings
R_ = [-18.3921, -15.1886, -12.2722, -9.6468, -7.3198, -5.2925, -3.5710, -2.164, -1.0753, -0.3050]
Prad_sim = [5.1086, 4.5802, 4.0441, 3.5318, 3.0742, 2.7724, 2.6058, 2.4659, 2.2395, 1.7829]
ne_sim = [9.1732, 8.7446, 8.3351, 7.8059, 6.9740, 6.2106, 5.5308, 4.9412, 4.4436, 4.0410] * 1e19
c_z_sim = [0.1944, 0.2034, 0.2127, 0.2271, 0.2583, 0.2992, 0.3524, 0.4204, 0.5046, 0.6109] *0.01
Te_sim = [1538.9561, 1359.4824, 1193.0386, 1038.3017, 891.58594, 750.65806, 617.58147, 496.69831, 393.1432, 312.39886]

Lz_Kr_sim = [Lz[:Kr].Lz_tot(ne_,Te_) for (ne_,Te_) in zip(ne_sim,Te_sim)]
plot(R_, Lz_Kr_sim)
plot(R_, Lz_Kr_sim .* c_z_sim .* ne_sim .^2 )



# plot(title="Average charge s