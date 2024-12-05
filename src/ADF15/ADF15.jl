#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

include("parser.jl")

Element = :W
Z = 0.0

lambda_input = 4008.8 # Angstrom
dens_input = vec([1e17 1e18 1e19 1e20]) # m⁻³
temp_input = vec([25.0 25.0 25.0 25.0]) # eV

datafile = get_adf15_datafile(Element, Z; bundling_model="ic")
log10pec_dict, meta = read_adf15(datafile)

PEC, lambda = get_interpolated_value(log10pec_dict, lambda_input, dens_input, temp_input)

pec_struct = get_photon_emissivity_coeff(datafile)