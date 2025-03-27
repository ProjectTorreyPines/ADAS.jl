#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#


# @__DIR__:
# This macro represents the directory of the file in which this code resides.
# For example, if this file is located at /home/user/MyModule/src/myfile.jl, @__DIR__ evaluates to /home/user/MyModule/src.ADASPaths = Dict{Symbol,Dict{Symbol,String}}
const data_paths = Dict(:parsed_data=>Dict(), :raw_data=>Dict())


"""
    set_adas_raw_data_directory!(dir=joinpath(@__DIR__, "../"))

Set the root directory for ADAS raw data files (adf11 and adf15).

The function iterates through the ADAS data types (`:adf11` and `:adf15`),
and sets the corresponding paths in the `data_paths[:raw_data]` dictionary.
If the specified directory does not exist, it creates the directory.

# Arguments
- `dir::String`: The root directory where the "data" subdirectory containing "adf11" and "adf15" subdirectories are located.
                   Defaults to the parent directory of the current file.
"""
function set_adas_raw_data_directory!(dir=joinpath(@__DIR__, "../"))
    for adas_type in [:adf11,:adf15]
        data_paths[:raw_data][adas_type] = joinpath(dir, "data/$(string(adas_type))")
    if !isdir(data_paths[:raw_data][adas_type])
            mkpath(data_paths[:raw_data][adas_type])
    end
end
end


"""
    set_adas_parsed_data_directory!(dir::String=joinpath(@__DIR__, "../"))

Set the directory where parsed ADAS data files are stored.

This function updates the `data_paths` dictionary with the specified directory for both `adf11` and `adf15` data types.
If the specified directory does not exist, it creates the directory.

# Arguments
- `dir::String`: The base directory where the parsed data directories for `adf11` and `adf15` are located.
                   Defaults to the parent directory of the current file.
"""
function set_adas_parsed_data_directory!(dir=joinpath(@__DIR__, "../"))
    for adas_type in [:adf11, :adf15]
        data_paths[:parsed_data][adas_type] = joinpath(dir, "parsed_data/$(string(adas_type))")
        if !isdir(data_paths[:raw_data][adas_type])
            mkpath(data_paths[:raw_data][adas_type])
        end
    end
end

set_adas_raw_data_directory!()
set_adas_parsed_data_directory!()

mutable struct ADASDatabase
    adf11::Dict
    adf15::Dict
end
Base.getindex(db::ADAS.ADASDatabase, s::Symbol) = getproperty(db, s)
const ADASdb = ADASDatabase(Dict(),Dict())

(db::ADASDatabase)(element, adas_type) = retrieve_element_database(element, db, adas_type)

function retrieve_element_database(element::Union{String,Symbol}, db::ADASDatabase, adas_type::Union{Symbol,String})
    element = string(element)
    element = lowercase(element)
    adas_type = Symbol(adas_type)
    @assert (adas_type ∈ propertynames(db)) "getting data of type $adas_type is not implemented yet... Available types: $(propertynames(db))"
    data = getproperty(db, adas_type)
    @debug "looking for $element in $(keys(data))"

    if element ∉ keys(data)
        load_element_data!(data, element, adas_type)
    end

    @assert element ∈ keys(data) "Cannot find element '$element' in database. Elements available are: $(collect(keys(data)))"
    return data[element]
end

function retrieve_element_data(data_elem::Dict; adas_type=missing, kwargs...)
    if adas_type == :adf11
        return retrieve_adf11_element_data(data_elem; kwargs...)
    elseif adas_type == :adf15
        return retrieve_adf15_element_data(data_elem; kwargs...)
    else
        error("Getting data of type $adas_type is not implemented yet...")
    end
end

"""alias for retrieve_element_data"""
retrieve_ADAS_data(args...; kw...) = retrieve_element_data(args...; kw...)

function retrieve_element_data(element::Union{Symbol,String}; adas_type=missing, type="scd", kw...)
    if ismissing(adas_type)
        adas_type = get_adas_type(type) 
    end
    @assert istype(adas_type,type) "Cannot find type $type for adas_type $adas_type. Available types are: $(collect(keys(ADASdb[adas_type])))"
    retrieve_element_data(ADASdb(element, adas_type); adas_type, type, kw...)
end 

function load_element_data!(data, element, adas_type)
    data[element] = get_element_data(element, adas_type)
end
function get_element_data(element::String, adas_type::Symbol)::Dict
    filepath = get_parsed_data_filepath(element, data_paths[:parsed_data][adas_type])
    if !isfile(filepath)
        make_database(element, paths, adas_type)
    end
    println("Getting element $element data from file: $filepath")
    data = load_data(filepath)
    return data
end

function get_adas_type(type)
    String(type) ∈ keys(adas_type_dict[:adf15]) && return :adf15
    String(type) ∈ keys(adas_type_dict[:adf11]) && return :adf11
    error("Cannot find type $type in available types: \n adf11: $(collect(keys(adas_type_dict[:adf11]))) \n adf15: $(collect(keys(adas_type_dict[:adf15])))")
end

istype(adas_type, type) = String(type) ∈ collect(keys(adas_type_dict[adas_type]))





export show_ADAS_files



function make_database(adas_type::Symbol)
    println("building database for $adas_type files in $(data_paths[:parsed_data][adas_type])")
    data = get_file_list(data_paths[:raw_data][adas_type], adas_type)
    dict = vec2dict(data)
    return dump_data(data_paths[:parsed_data][adas_type], dict)
end
function make_database(adas_types::Vector{Symbol})
    for adas_type in adas_types
        make_database(adas_type)
    end 
end

function make_database(element, adas_type::Symbol)
    element = lowercase(element)
    data = get_file_list(data_paths[:raw_data][adas_type], adas_type)
    @assert element ∈ keys(data)
    return dump_data(element, data_paths[:parsed_data][adas_type], data)
end

build_ADAS_database(; adas_type=[:adf15, :adf11]) = make_database(adas_type)


# function build_ADAS_database(paths, adas_type)

#     data = get_database(paths[:raw_data][adas_type], adas_type)
#     dump_data(paths[:parsed_data][adas_type], data)
# end


