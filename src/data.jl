#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#

ADASPaths = Dict{Symbol,Dict{Symbol,String}}
const parsed_data_directory = Dict(:adf11 => joinpath(@__DIR__, "../parsed_data/adf11"))

const data_directory = Dict(:adf11 => joinpath(@__DIR__, "../data/adf11"))

mutable struct ADASData
    adf11::ADF11
    paths::ADASPaths
end

const ADASdata = ADASData(Dict(), Dict(:parsed_data => parsed_data_directory, :raw_data => data_directory))

(db::ADASData)(element; adas_type::Symbol=:adf11) = (db::ADASData)(element, adas_type)

(db::ADASData)(element, adas_type::String) = (db::ADASData)(element, Symbol(adas_type))

function (ADASdata::ADASData)(element, adas_type::Symbol)
    @assert adas_type == :adf11 "getting data of type $adas_type is not implemented yet..."
    return retrieve_element_data(element, getfield(ADASdata, adas_type), adas_type)
end

retrieve_element_data(element::Symbol, data, adas_type::Symbol) = retrieve_element_data(string(element), data, adas_type)

function retrieve_element_data(element::String, data, adas_type::Symbol)
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

function get_element_data(element::String, paths::ADASPaths, adas_type::Symbol)
    filepath = get_data_filepath(element, paths[:parsed_data][adas_type])
    println("Getting element $element data from file: $filepath")
    if !isfile(filepath)
        make_database(element, paths, adas_type)
    end
    data = load_data(filepath)::Dict
    return data
end

const colors = Dict{Int,Symbol}()
colors[0] = :blue
colors[1] = :red
colors[2] = :green
colors[3] = :magenta

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

get_data_filepath(element, directory) = joinpath(directory, "$element.jld2")

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

get_database(directory::String, adas_type) = get_database(directory::String, ADASType{adas_type}())

get_database(directory::String, ::ADASType{:adf11}) = get_database(get_list_file(directory, adf11File))

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