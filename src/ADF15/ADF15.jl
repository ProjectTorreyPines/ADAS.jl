#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

include("parser.jl")

path = "/home/cappellil/adf15_python/pec40#w_cl#w1.dat"
order = Int(1)

lambda_input = 142.0 # Angstrom
dens_input = 1e15 # cm⁻³
temp_input = 20.0 # eV

PEC = get_interpolated_value(log10pec_dict, lambda_input, dens_input, temp_input)