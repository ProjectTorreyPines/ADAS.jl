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

adas_type_dict[:adf11] = Dict(
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
    for (k, v) in adas_type_dict[:adf11]
        println("'$k' : $(v[2])")
    end
    println("--------------------------------")
end

function get_type_adf11(type)
    @assert type ∈ keys(adas_type_dict[:adf11]) "$type not in $(keys(adas_type_dict[:adf11]))"
    return adas_type_dict[:adf11][type][1]
end

struct ADF11Block <: ADASBlock
    header::String
    content::Vector{String}
end
ADF11Block(s::Vector{String}) = ADF11Block(s[1], s[2:end])
struct ADF11Header <: ADASHeader
    n_ions::Int64
    n_ne::Int64
    n_Te::Int64
    imin::Int64
    imax::Int64
    details::String
end

function ADF11Header(block::ADF11Block)
    header = block.content[1]
    data = parse.(Int64, split(header)[1:5])
    details = join(split(header)[6:end], " ")
    @debug "ADF11Header : data=$data; details=$details"
    return ADF11Header(data..., details)
end

struct ADASAxis
    log10_Te::Vector{Float64}
    log10_ne::Vector{Float64}
    Te::Vector{Float64}
    ne::Vector{Float64}
    unit_Te::String
    unit_ne::String
end

function ADASAxis(block::ADF11Block, header::ADF11Header; unit_ne="m^-3") #TODO: use unitful to manage units
    log10_Te = Vector{Float64}()
    log10_ne = Vector{Float64}()
    idx = 1
    while length(log10_ne) < header.n_ne
        push!(log10_ne, parse.(Float64, [n for n in split(block.content[idx])])...)
        idx += 1
    end
    while length(log10_Te) < header.n_Te
        push!(log10_Te, parse.(Float64, [t for t in split(block.content[idx])])...)
        idx += 1
    end
    if unit_ne == "m^-3"
        log10_ne .= log10_ne .+ 6.0
    elseif unit_ne == "cm^-3"

    else
        error("Invalid unit_ne: $unit_ne")
    end
    return ADASAxis(log10_Te, log10_ne, 10 .^ log10_Te, 10 .^ log10_ne, "eV", unit_ne)
end

struct ADF11Rates{T} <: ADASRates{T}
    igrnd::Union{Int64,Nothing}
    iptr::Union{Int64,Nothing}
    Z::Int64
    log10_values::T
    values::T
end

function ADF11Rates(block::ADF11Block, header::ADF11Header; unit_rates="m^3")
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
    return ADF11Rates(igrd, iptr, Z, log10_values, 10.0 .^ log10_values)
end

struct adf11Data{T}
    filepath::String
    header::ADF11Header
    axis::ADASAxis
    rates::T
    units::String
end

function adf11Data(filepath::String; metastable=false)
    lines = read_file(filepath)
    comments = get_comments(lines)
    units = get_units(comments)
    blocks = split_blocks(lines)
    header = ADF11Header(blocks[1])
    if metastable
        axis = ADASAxis(blocks[3], header)
        rates = Dict{Int64,Vector{ADF11Rates}}()
        for block in blocks
            Z = get_block_attr(block.header, "Z1")
            if Z !== nothing
                if !(Z in collect(keys(rates)))
                    rates[Z] = Vector{ADF11Rates}()
                end
                push!(rates[Z], ADF11Rates(block, header))
            end
        end
    else
        axis = ADASAxis(blocks[2], header)
        rates = Dict{Int64,ADF11Rates}()
        for block in blocks
            Z = get_block_attr(block.header, "Z1")
            if Z !== nothing
                rates[Z] = ADF11Rates(block, header)
            end
        end
    end
    return adf11Data(filepath, header, axis, rates, units)
end

struct adf11File{T} <: ADASFile{T}
    element::String
    filepath::String
    year::String
    metastable::Bool
    type::String
    data::adf11Data
    md5::Vector{UInt8}
end

function adf11File(filepath::String)
    filename = basename(filepath)
    element = match(Regex("(?<=\\_)(.*?)(?=\\.|\\#)"), filename).captures[1]
    type = match(Regex("(.*?)(?=[0-9]{2})"), filename).captures[1]
    year = match(Regex("([0-9]{2})"), filename).captures[1]
    if match(Regex("(?<=[0-9]{2})(r)(?=\\_)"), filename) !== nothing
        metastable = true
    else
        metastable = false
    end
    data = adf11Data(filepath; metastable=metastable)
    return adf11File{get_type_adf11(type)}(element, filepath, year, metastable, type, data, checksum(filepath))
end



const ADF11 = Dict{String,Dict{String,Dict{String,Dict{String,adf11File}}}}