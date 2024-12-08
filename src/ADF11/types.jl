#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#


abstract type Ionization <: ADASRate end
abstract type Recombination <: ADASRate end
abstract type ContinuumRadiation <: ADASRate end
abstract type LineRadiation <: ADASRate end
abstract type ChargeExchange <: ADASRate end
abstract type ChargeExchangeRadiation <: ADASRate end
abstract type SXRLineRadiation <: ADASRate end
abstract type SXRContinuumRadiation <: ADASRate end
abstract type RecombinationBremsstrahlung <: ADASRate end
abstract type SXRSensitivity <: ADASRate end
abstract type SXRBremsstrahlung <: ADASRate end
abstract type Bremsstrahlung <: ADASRate end
abstract type MeanIonisationPotential <: ADASRate end
abstract type CrossCouplingCoeffs <: ADASRate end
abstract type ParentCrossCouplingCoeffs <: ADASRate end
abstract type MeanChargeStateSquared <: ADASRate end
abstract type MeanChargeState <: ADASRate end

const adf11_types = Dict(
    "acd" => (Recombination, "effective recombination coefficients"),
    "scd" => (Ionization, "effective ionization coefficients"),
    "prb" => (RecombinationBremsstrahlung, "Continuum and line power driven by recombination and Bremsstrahlung of dominant ions"),
    "plt" => (LineRadiation, "Line power driven by excitation of dominant ions"),
    "ccd" => (ChargeExchange, "Charge exchange effective recombination coefficients (with D)"),
    "prc" => (ChargeExchangeRadiation, "Line power due to charge transfer from thermal neutral hydrogen to dominant ions (Charge exchange emission)"),
    "pls" => (SXRLineRadiation, "Line power from selected transitions of dominant ions"),
    "prs" => (SXRContinuumRadiation, "continuum radiation in the SXR range"),
    "brs" => (Bremsstrahlung, "continuum spectral bremstrahlung"),
    "fis" => (SXRSensitivity, "sensitivity in the SXR range"),
    "pbs" => (SXRBremsstrahlung, "impurity bremsstrahlung in SXR range, also included in prs files"),
    "ecd" => (MeanIonisationPotential, "Mean Ionisation Potential"),
    "qcd" => (CrossCouplingCoeffs, "Cross-coupling coefficients"),
    "xcd" => (ParentCrossCouplingCoeffs, "Parent Cross-coupling coefficients"),
    "ycd" => (MeanChargeStateSquared, "Mean Charge State Squared"),
    "zcd" => (MeanChargeState, "Mean Charge State")
)

function show_adf11_types()
    println("----------- adf11 types ------------")
    for (k, v) in adf11_types
        println("'$k' : $(v[2])")
    end
    return print("--------------------------------------")
end

function get_type_adf11(type)
    @assert type ∈ keys(adf11_types) "$type not in $(keys(adf11_types))"
    return adf11_types[type][1]
end

struct ADASBlock
    header::String
    content::Vector{String}
end

struct ADASHeader
    n_ions::Int64
    n_ne::Int64
    n_Te::Int64
    imin::Int64
    imax::Int64
    details::String
end

function ADASHeader(block::ADASBlock)
    header = block.content[1]
    data = parse.(Int64, split(header)[1:5])
    details = join(split(header)[6:end], " ")
    @debug "ADASHeader : data=$data; details=$details"
    return ADASHeader(data..., details)
end

struct ADASAxis
    log_Te::Vector{Float64}
    log_ne::Vector{Float64}
    Te::Vector{Float64}
    ne::Vector{Float64}
    unit_Te::String
    unit_ne::String
end

function ADASAxis(block::ADASBlock, header::ADASHeader; unit_ne="m^-3") #TODO: use unitful to manage units
    log_Te = Vector{Float64}()
    log_ne = Vector{Float64}()
    idx = 1
    while length(log_ne) < header.n_ne
        push!(log_ne, parse.(Float64, [n for n in split(block.content[idx])])...)
        idx += 1
    end
    while length(log_Te) < header.n_Te
        push!(log_Te, parse.(Float64, [t for t in split(block.content[idx])])...)
        idx += 1
    end
    if unit_ne == "m^-3"
        log_ne .= log_ne .+ 6.0
    end
    return ADASAxis(log_Te, log_ne, 10 .^ log_Te, 10 .^ log_ne, "eV", unit_ne)
end

struct ADASRates{T}
    igrnd::Union{Int64,Nothing}
    iptr::Union{Int64,Nothing}
    Z::Int64
    log_values::T
    values::T
end

function ADASRates(block::ADASBlock, header::ADASHeader; unit_rates="m^3")
    Z = get_block_attr(block.header, "Z1")
    igrd = get_block_attr(block.header, "IGRD")
    iptr = get_block_attr(block.header, "IPRT")
    tmp = Vector{Float64}() #zeros(Float64,header.n_ne * header.n_Te)
    idx = 1
    while length(tmp) < header.n_ne * header.n_Te
        push!(tmp, parse.(Float64, [replace(r, "D" => "e") for r in split(block.content[idx])])...)
        idx += 1
    end
    log10_values = reshape(tmp, header.n_ne, header.n_Te)
    @debug "igrd = $igrd; iptr = $iptr; Z = $Z"
    if unit_rates == "m^3"
        log10_values .= log10_values .- 6.0
    end
    return ADASRates(igrd, iptr, Z, log10_values, 10.0 .^ log10_values)
end

struct adf11Data{T}
    filepath::String
    header::ADASHeader
    axis::ADASAxis
    rates::T
    units::String
end

function adf11Data(filepath::String; metastable=false)
    lines = read_adas_file(filepath)
    comments = get_comments(lines)
    units = get_units(comments)
    blocks = split_blocks(lines)
    header = ADASHeader(blocks[1])
    if metastable
        axis = ADASAxis(blocks[3], header)
        rates = Dict{Int64,Vector{ADASRates}}()
        for block in blocks
            Z = get_block_attr(block.header, "Z1")
            if Z !== nothing
                if !(Z in collect(keys(rates)))
                    rates[Z] = Vector{ADASRates}()
                end
                push!(rates[Z], ADASRates(block, header))
            end
        end
    else
        axis = ADASAxis(blocks[2], header)
        rates = Dict{Int64,ADASRates}()
        for block in blocks
            Z = get_block_attr(block.header, "Z1")
            if Z !== nothing
                rates[Z] = ADASRates(block, header)
            end
        end
    end
    return adf11Data(filepath, header, axis, rates, units)
end

struct adf11File{T} <: ADASFile{T}
    name::String
    element::String
    path::String
    year::String
    metastable::Bool
    type::String
    data::adf11Data
    md5::Vector{UInt8}
end

function adf11File(filename::String, filepath::String)
    name = filename
    element = match(Regex("(?<=\\_)(.*?)(?=\\.|\\#)"), filename).captures[1]
    type = match(Regex("(.*?)(?=[0-9]{2})"), filename).captures[1]
    year = match(Regex("([0-9]{2})"), filename).captures[1]
    if match(Regex("(?<=[0-9]{2})(r)(?=\\_)"), filename) !== nothing
        metastable = true
    else
        metastable = false
    end
    data = adf11Data(filepath; metastable=metastable)
    return adf11File{get_type_adf11(type)}(name, element, filepath, year, metastable, type, data, checksum(filepath))
end

function checksum(filepath)
    open(filepath) do f # The open function opens the file located at filepath and assigns the file handle to the variable f within the scope of the do block.
                        # When the block finishes executing, open automatically closes the file, even if an error occurs, ensuring safe and efficient file handling.
        return Array(md5(f)) # Converts the hash into an array of bytes for easy manipulation or comparison.
                             # MD5 (Message Digest 5) is one specific hash function. It produces a fixed-size output (e.g., 128 bits for MD5), regardless of the size of the input.
                             # Is deterministic, meaning the same input will always produce the same hash.
                             # For MD5 the hash is represented with 32 hexadecimal characters (e.g., 5d41402abc4b2a76b9719d911017c592)
                             # each hash is a unique blue print uniquely identifying data
    end
end

const ADF11 = Dict{String,Dict{String,Dict{String,Dict{String,adf11File}}}}