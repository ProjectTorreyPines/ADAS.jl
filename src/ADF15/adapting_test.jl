
using Logging
using MD5
using AbstractTrees
using FileIO
using Downloads
using Interpolations


abstract type ADASFile{T} end
struct ADASType{T} end

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


const ADF11 = Dict{String,Dict{String,Dict{String,Dict{String,adf11File}}}}


## ----- ADF15 ----- ##

struct adf15Data{T}
    filepath::String
    header::ADASHeader
    axis::ADASAxis
    rates::T
    units::String
end

function adf15Data(filepath::String; metastable=false)
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

struct adf15File{T} <: ADASFile{T}
    name::String
    element::String
    path::String
    year::String
    metastable::Bool
    type::String
    data::adf11Data
    md5::Vector{UInt8}
end

function adf15File(filename::String, filepath::String)
    name = filename
    element = match(Regex("(?<=\\_)(.*?)(?=\\.|\\#)"), filename).captures[1]
    type = match(Regex("(.*?)(?=[0-9]{2})"), filename).captures[1]
    year = match(Regex("([0-9]{2})"), filename).captures[1]
    if match(Regex("(?<=[0-9]{2})(r)(?=\\_)"), filename) !== nothing
        metastable = true
    else
        metastable = false
    end
    data = adf15Data(filepath; metastable=metastable)
    return adf15File{get_type_adf11(type)}(name, element, filepath, year, metastable, type, data, checksum(filepath))
end

const ADF15 = Dict{String,Dict{String,Dict{String,Dict{String,adf15File}}}}

# paths for data given a priori

const parsed_data_directory = Dict(:adf11 => joinpath(@__DIR__, "../parsed_data/adf11"), :adf15 => joinpath(@__DIR__, "../parsed_data/adf15"))
const data_directory = Dict(:adf11 => joinpath(@__DIR__, "../data/adf11"), :adf15 => joinpath(@__DIR__, "../data/adf15"))

const ADASPaths = Dict{Symbol, Dict{Symbol, String}}

mutable struct ADASData
    adf11::ADF11
    adf15::ADF15 # adding ADF15 the methods won't work
    paths::ADASPaths
end

(db::ADASData)(element; adas_type::Symbol=:adf11) = (db::ADASData)(element, adas_type)

(db::ADASData)(element, adas_type::String) = (db::ADASData)(element, Symbol(adas_type))

function (ADASdata::ADASData)(element, adas_type::Symbol)
    @assert adas_type == :adf11  || adas_type == :adf15 "getting data of type $adas_type is not implemented yet..."
    return retrieve_element_data(element, getfield(ADASdata, adas_type), adas_type)
end

const ADASdata = ADASData(Dict(), Dict(), Dict(:parsed_data => parsed_data_directory, :raw_data => data_directory))

## ----- functions to retrieve data from paths ----- ## 

get_data_filepath(element, directory) = joinpath(directory, "$element.jld2")


function get_element_data(element::String, paths::ADASPaths, adas_type::Symbol)::Dict
    filepath = get_data_filepath(element, paths[:parsed_data][adas_type])
    println("Getting element $element data from file: $filepath")
    if !isfile(filepath)
        make_database(element, paths, adas_type)
    end
    data = load_data(filepath)
    return data
end


function make_database(paths, adas_type)
    println("building database for $adas_type files in $paths")
    data = get_database(paths[:raw_data][adas_type], adas_type)
    return dump_data(paths[:parsed_data][adas_type], data)
end

function make_database(element, paths, adas_type)
    element = lowercase(element)
    data = get_database(paths[:raw_data][adas_type], adas_type)
    @assert element ∈ keys(data)
    return dump_data(element, paths[:parsed_data][adas_type], data)
end

# ADASType{T} is an empty structure of type T

# you can use get_database either by defining directory path and adas_type or path and adf11File, it's the same
get_database(directory::String, adas_type::Symbol) = get_database(directory::String, ADASType{adas_type}())

get_database(directory::String, adf11File::ADASType{:adf11}) = get_database(get_list_file(directory, adf11File))

get_database(directory::String, adf11File::ADASType{:adf15}) = get_database(get_list_file(directory, adf11File))


function get_database(adf11_files::Vector{<:adf11File})
    database = Dict{String,Dict{String,Dict{String,Dict{String,adf11File}}}}()
    for file in adf11_files
        if !(file.element in keys(database))
            database[file.element] = Dict{String,Dict{String,Dict{String,adf11File}}}()
        end
        if !(file.type in keys(database[file.element]))
            database[file.element][file.type] = Dict{String,Dict{String,adf11File}}()
        end
        if !(file.year in keys(database[file.element][file.type]))
            database[file.element][file.type][file.year] = Dict{String,adf11File}()
        end
        if file.metastable
            database[file.element][file.type][file.year]["metastable"] = file
        else
            database[file.element][file.type][file.year]["ground"] = file
        end

    end
    return database
end


function get_list_file(directory, FileType)
    files_ = Vector{FileType}()
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if occursin(".dat", file)
                push!(files_, FileType(file, joinpath(root, file)))
            end
        end
    end
    return files_
end



retrieve_element_data(element::Symbol, data, adas_type::Symbol) = retrieve_element_data(string(element), data, adas_type)

function retrieve_element_data(element::String, data::adf11File, adas_type::Symbol)
    element = string(element)
    element = lowercase(element)
    @debug "looking for $element in $(keys(data))"

    if element ∉ keys(data)
        println("ADASdata.paths:", ADASdata.paths)
        data[element] = get_element_data(element, ADASdata.paths, adas_type)[element]
    end

    @assert element ∈ keys(data) "Cannot find element '$element' in database. Elements available are: $(collect(keys(data)))"
    return data[element]
end

function retrieve_element_data(args...; adas_type::Symbol=:adf11, kwargs...)
    if adas_type == :adf11
        return retrieve_adf11_element_data(args...; kwargs...)
    else
        error("Getting data of type $adas_type is not implemented yet...")
    end
end

function retrieve_ADAS_data(element::Symbol; kw...)
    return retrieve_ADAS_data(string(element); kw...)
end

function retrieve_ADAS_data(element::String; year::Union{String,Missing}="latest", type::String="scd", metastable::Bool=false, adas_type=:adf11)
    return retrieve_element_data(ADASdata(element, adas_type); year=year, type=type, metastable=metastable, adas_type=adas_type)
end



const colors = Dict{Int,Symbol}()
colors[0] = :blue
colors[1] = :red
colors[2] = :green
colors[3] = :magenta

# The code defines how to display ADAS data in a tree-like structure using the AbstractTrees package

show_ADAS_data(element; adas_type=:adf11, kw...) = AbstractTrees.print_tree(ADASdata(element, adas_type); kw...)
AbstractTrees.printnode(io::IO, a::ADASData) = printstyled("ADAS data: $(a.adf11.element)"; bold=true)
AbstractTrees.printnode(io::IO, a::Dict{String,Dict{String,adf11File}}) = printstyled("years")
AbstractTrees.printnode(io::IO, a::Dict{String,Dict{String,Dict{String,ADAS.adf11File}}}) = printstyled("adf11")
AbstractTrees.printnode(io::IO, a::Dict{String,ADAS.adf11File}) = nothing


show_ADAS_data(; adas_type=:adf11, kw...) = AbstractTrees.print_tree(get_database(ADASdata.paths[:raw_data][adas_type], adas_type); kw...)


function dump_data(element, directory, data)
    file_path = get_data_filepath(element, directory)
    FileIO.save(file_path, data)
    return nothing
end

function dump_data(directory, data)
    for (el, d) in data
        dump_data(el, directory, d)
    end
    return nothing
end



function load_data(path::String)
    @debug "Loading jld2 file $path"
    return FileIO.load(path) #BSON.load(path, @__MODULE__)
end

load_data(element::String, path::String) = load_data(get_data_filepath(element, path))

function read_adas_file(filepath::String)
    @debug "Reading ADAS file: $filepath"
    @assert isfile(filepath) "File not found: $filepath"
    fid = open(filepath, "r")
    lines = readlines(fid)
    close(fid)
    return lines
end



function build_ADAS_database(; parsed_data_path::Union{Missing,String}=missing)
    local_ADASdata = deepcopy(ADASdata)
    if !ismissing(parsed_data_path)
        for (k, v) in local_ADASdata.paths[:parsed_data]
            local_ADASdata.paths[:parsed_data][k] = joinpath(parsed_data_path, string(k))
        end

    end

    for adas_type in (f for f in propertynames(local_ADASdata) if f != :paths)
        make_database(local_ADASdata.paths, adas_type)
    end
end

# function build_ADAS_database(paths, adas_type)

#     data = get_database(paths[:raw_data][adas_type], adas_type)
#     dump_data(paths[:parsed_data][adas_type], data)
# end


function Base.show(io::IO, ::MIME"text/plain", data::ADASData)
    return print(io, "ADAS data: adf11: $(collect(keys(data.adf11)))")
end

function Base.show(io::IO, data::ADASData)
    return print(io, "ADAS data: adf11: $(join(collect(keys(data.adf11)),"; "))")
end

Base.show(io::IO, ::MIME"text/plain", file::adf11File) = AbstractTrees.print_tree(file; maxdepth=1, indicate_truncation=false)
AbstractTrees.printnode(file::adf11File) = printstyled(io, "ADAS adf11 data: $(file.name)"; bold=true)
AbstractTrees.children(file::adf11File) = Dict(f => getproperty(file, f) for f in propertynames(file) if f != :data && f != :md5)


function Base.show(io::IO, data::adf11File)
    return print(io, "ADAS data | adf11 : $(data.name) ")
end