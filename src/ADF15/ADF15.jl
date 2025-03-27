#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#
include("types.jl")



function vec2dict(adf15_files::Vector{<:adf15File})
    database = Dict()
    for file in adf15_files
        if !(file.element in keys(database))
            database[file.element] = Dict()
        end
        if !(file.Z in keys(database[file.element]))
            database[file.element][file.Z] = Dict()
        end
        if !(file.year in keys(database[file.element][file.Z]))
            database[file.element][file.Z][file.year] = Dict()
        end
        if !(file.model in keys(database[file.element][file.Z][file.year]))
            database[file.element][file.Z][file.year][file.model] = file
        end

    end
    return database
end


function retrieve_adf15_element_data(data::Dict; Z=missing, year::Union{String,Missing}="latest", type="pec", model=missing)
    @assert !(Z isa Missing) "Z must be provided to retrieve adf15 data ... Use keyword Z=XX where XX is the atomic number"
    Z = "$Z"
    @assert Z ∈ keys(data) "Cannot find element '$Z' in database. Charge state available are: $(collect(keys(data)))"
    data[Z]
    if (year isa Missing) || (year == "latest")
            year = sort!(collect(keys(data[Z])))[end]
    end
    @assert year ∈ keys(data[Z]) "Cannot find year '$year' in data. Years avaible are: $(collect(keys(data[Z])))"

    @assert !(model isa Missing) "Model must be provided to retrieve adf15 data ... Use keyword model=XX where XX is the model. 
\n Available models for Z=$Z are: $(collect(keys(data[Z][year])))"
    model = lowercase(model)
    @assert model ∈ keys(data[Z][year]) "Cannot find model '$model' in database. Models available are: $(collect(keys(data[Z][year])))"
  
    return data[Z][year][model]
end


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

function parse_info(block)
    lines = split(block, '\n')
    parsed_data = []

    for line in lines
        # Skip comment lines or empty lines
        if startswith(line, "C") || isempty(strip(line))
            continue
        end

        # Extract fields using a regular expression
        match = match(r"^\s*(\d+)\.\s+([\d.]+)\s+([\d()\-.\s]+)\s+(\w+)\s+([\w\s#]*)\s+([T\s]*)\s+(\d*)$", line)
        if match !== nothing
            push!(parsed_data, (
                ISEL = parse(Int, match.captures[1]),
                WAVELENGTH = parse(Float64, match.captures[2]),
                TRANSITION = strip(match.captures[3]),
                TYPE = strip(match.captures[4]),
                METASTABLE = strip(match.captures[5]),
                IMET = strip(match.captures[6]),
                NMET = isempty(match.captures[7]) ? missing : parse(Int, match.captures[7])
            ))
        end
    end

    return parsed_data
end