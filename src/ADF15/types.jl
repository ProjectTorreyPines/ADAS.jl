#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
Contributor : Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

abstract type Excitation <: ADASRate end


abstract type PhotonEmissivityCoeff <: ADASRate end

adas_type_dict[:adf15] = Dict(
    "pec" => (PhotonEmissivityCoeff, "Emission related to excitation and relaxation"),
    # "rec" => (Recombination, "Free electron binds to ion, emitting photon"),
    # "cx" => (ChargeExchange, "Emission caused by a neutral atom swapping an electron with an ion") 
    )


struct ADF15Block
    header::Dict
    content::Vector{String}
end
ADF15Block(s::Vector{String}) = ADF15Block(adf15_block_header(s[1]), s[2:end])


struct ADF15Header <: ADASHeader
    wavelengths::Vector{Float64}
    n_ne::Int64
    n_Te::Int64
    info::Dict
end

function ADF15Header(blocks::Vector{ADF15Block})
    @assert all((b.header["dens_dim"] for b in blocks) .== blocks[1].header["dens_dim"])
    @assert all((b.header["temp_dim"] for b in blocks) .== blocks[1].header["temp_dim"])
    n_Te = blocks[1].header["temp_dim"]
    n_ne = blocks[1].header["dens_dim"]
    wavelengths = [b.header["wavelength"] for b in blocks]
    info = Dict(b.header["wavelength"] => b for b in blocks)
    header = ADF15Header(wavelengths, n_ne, n_Te, info)
    return header
end

struct ADF15Rates{T} <: ADASRates{T}
    wavelength::Float64
    log10_values::T
    values::T
    info::Dict
    unit_rates::String
end


function ADF15Rates(block::ADF15Block, header::ADF15Header; unit_rates="m^3.s^-1")
    wavelength = block.header["wavelength"]
    info = block.header
    tmp = Vector{Float64}() #zeros(Float64,header.n_ne * header.n_Te)
    idx = 1
    s = 0
    while s < header.n_ne + header.n_Te
        s = s + length(split(block.content[idx]))
        idx += 1
    end

    while length(tmp) < header.n_ne * header.n_Te
        push!(tmp, parse.(Float64, [replace(r, "D" => "e") for r in split(block.content[idx])])...)
        idx += 1
    end
    values = Matrix(transpose(reshape(tmp, header.n_Te, header.n_ne)))
    if unit_rates == "m^3.s^-1"
        values .= values * 1e-6
    elseif unit_rates == "cm^3.s^-1"
    else
        error("Invalid unit_rates: $unit_rates")
    end
    return ADF15Rates(wavelength, log10.(values), values, info, unit_rates)
end

struct adf15Data{T}
    filepath::String
    header::ADF15Header
    axis::ADASAxis
    rates::T
    units::String
end

function adf15Data(fp::String; unit_rates="m^3.s^-1", unit_ne="m^-3")
    lines = ADAS.read_file(fp)
    comments = ADAS.get_comments(lines)
    data = filter(line -> !(startswith(line, "C") || startswith(line, "c") || startswith(line, "//") || line == "" || all(isspace, line)), lines)
    nblocks = parse(Int, match(r"\d+", data[1]).match)
    last_index = findlast(line -> startswith(line, "---"), data)
    if last_index isa Nothing
        last_index = 1
    end 
    d = data[last_index+1:end]
    lblock = Int64(length(d) / nblocks)
    @assert length(data[last_index+1:end]) % nblocks == 0
    blocks = [ADF15Block(d[i*lblock+1:(i+1)*lblock]) for i in 0:nblocks-1]
    header = ADAS.ADF15Header(blocks)
    axis = ADAS.ADASAxis(blocks[1], header; unit_ne)
    rates = Dict(block.header["wavelength"] => ADF15Rates(block, header; unit_rates) for block in blocks)
    return adf15Data(fp, header, axis, rates, unit_rates)
end


function adf15_block_header(header::String)
    info = Dict()
    split_parts = split(header, r"[/:]")
    split_parts[1] = replace(split_parts[1], r"[a-zA-Z]" => "")
    numbers = filter(x -> !isempty(x), split(split_parts[1]))
    info["wavelength"] = parse(Float64, numbers[1])
    info["dens_dim"] = parse(Int, numbers[2])
    info["temp_dim"] = parse(Int, numbers[3])
    for s in split_parts[2:end]
        key, value = split(s, "="; limit=2)
        info[strip(key)] = strip(value)
    end
    return info
end


struct adf15File <: ADASFile{PhotonEmissivityCoeff}
    filepath::String
    element::String
    year::String
    model::String
    Z::String
    data::adf15Data
    md5::Vector{UInt8}
end


function ADASAxis(block::ADF15Block, header::ADF15Header; unit_ne="m^-3")
    Te = Vector{Float64}()
    ne = Vector{Float64}()
    idx = 1
    while length(ne) < header.n_ne
        push!(ne, parse.(Float64, [n for n in split(block.content[idx])])...)
        idx += 1
    end
    while length(Te) < header.n_Te
        push!(Te, parse.(Float64, [t for t in split(block.content[idx])])...)
        idx += 1
    end
    if unit_ne == "m^-3"
        ne .= ne .* 1e6
    end
    return ADASAxis(log10.(Te), log10.(ne), Te, ne, "eV", unit_ne)
end

"""
    get_adf15_filename_info(fp::String) -> Dict

Extracts information from the filename of an ADF15 file and returns it as a dictionary.

# Arguments
- `fp::String`: The file path of the ADF15 file.

# Returns
- `Dict`: A dictionary containing the following keys:
  - `"element"`: The chemical element extracted from the filename (e.g., "H", "He").
  - `"Z"`: The atomic number extracted from the filename.
  - `"model"`: The model name or identifier extracted from the filename.
  - `"year"`: The year associated with the PEC data extracted from the filename.

# Notes
The function assumes that the filename follows a specific naming convention, such as:
`#<element>_<model>#pec<year>_<Z>.dat`. If the filename does not match this pattern, the function may throw an error.
"""
function get_adf15_filename_info(fp::String)::Dict
    fn = basename(fp)
    info = Dict()
    info["element"] = match(r"#([a-zA-Z]+)_", fn).captures[end]
    info["Z"] = match(r"(\d+)\.dat", fn).captures[end]
    info["model"] = match(r"_(.*?)#", fn).captures[end]
    info["year"] = match(r"pec(\d+)", fn).captures[end]
    return info
end

function adf15File( fp::String)
    info = get_adf15_filename_info(fp)
    data = adf15Data(fp)
    return adf15File(fp, info["element"],  info["year"], info["model"], info["Z"], data, checksum(fp))
end


const ADF15 = Dict{String,Dict{String,Dict{String,Dict{String,adf15File}}}}
