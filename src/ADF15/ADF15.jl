#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

include("parser2.jl")

#= INPUT example and usage
Element = :W
Z = 0.0

lambda_input = 4008.8 # Angstrom
dens_input = vec([1e17 1e18 1e19 1e20]) # m⁻³
temp_input = vec([25.0 25.0 25.0 25.0]) # eV

PEC, lambda, pec_struct = get_pec(Element, Z)
=#

## ----- functions and output ----- ##

function get_pec(Element::Symbol, Z::Float64, lambda_input::Float64, dens_input::Vector{Float64}, temp_input::Vector{Float64}; bundling_model::String="ic")
    
    datafile = get_adf15_datafile(Element, Z; bundling_model)
    
    log10pec_dict, _ = read_adf15(datafile) 
    
    pec_struct = get_photon_emissivity_coeff(datafile)
    
    PEC, wavelength = get_pec_interpolated_value(log10pec_dict, lambda_input, dens_input, temp_input)
    
    return PEC, wavelength, pec_struct
end

# Dispatch for Element::String
get_pec(Element::String, kwargs...) = get_pec(Symbol(Element), kwargs...)