#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

include("parser.jl")

Element = :W
Z = 0.0

lambda_input = 4008.8 # Angstrom
dens_input = 1e20 # m⁻³
temp_input = 20.0 # eV

datafile = get_adf15_datafile(Element, Z)
log10pec_dict = read_adf15(datafile)
PEC = get_interpolated_value(log10pec_dict, lambda_input, dens_input, temp_input)