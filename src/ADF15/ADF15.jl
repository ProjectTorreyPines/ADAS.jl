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

struct PEC_struct{T,F,D}
    PEC::T
    wavelength::F
    pec_data::D
end

function get_pec(Element::Symbol, Z::Float64, lambda_input::Float64, dens_input::Union{Float64, Vector{Float64}}, temp_input::Union{Float64, Vector{Float64}}; bundling_model::String="ic")
    
    datafile = get_adf15_datafile(Element, Z; bundling_model)
    
    log10pec_dict, _ = read_adf15(datafile) 
    
    pec_struct = get_photon_emissivity_coeff(datafile)

    # if temp is a column vector transpose it
    if size(temp_input, 1) > size(temp_input, 2)
        temp = temp_input'  
    else
        temp = temp_input
    end

    # if dens is a row vector transpose it
    if size(dens_input, 1) < size(dens_input, 2)
        dens = dens_input'
    else 
        dens = dens_input  
    end
    
    PEC, wavelength = get_pec_interpolated_value(log10pec_dict, lambda_input, dens, temp)
    
    return PEC_struct(PEC, wavelength, pec_struct)
end

# Dispatch for Element::String
get_pec(Element::String, kwargs...) = get_pec(Symbol(Element), kwargs...)