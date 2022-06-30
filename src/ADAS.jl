module ADAS
using BSON
using RedefStructs
const default_adf11_database_file = joinpath(@__DIR__,"../data/adf11.bson")
const default_adf11_directory = joinpath(@__DIR__,"../data/adf11")
verbose = false #for debugging purpose
include("ADF11.jl")

mutable struct ADASDatabase
    adf11 :: Dict{Any,Any}
    function ADASDatabase()
        new(Dict())
    end
end

const ADAS_database = ADASDatabase()

# loading the database if it already exists otherwise build and save into bson
#adf11_database = get_adf11_database()

function retrieve_database(adas_type)
    if lowercase(adas_type) == "adf11"
        if length(ADAS_database.adf11) == 0
            ADAS_database.adf11 = get_adf11_database()
        end
        return ADAS_database.adf11
    else
        error("getting data of type $adas_type is not implemented yet...")
end 
end

function retrieve_data(args...; adas_type = "adf11", kwargs...)
    if adas_type == "adf11"
        return retrieve_adf11_data(args...;kwargs...)
    else
        error("getting data of type $adas_type is not implemented yet...")
    end
end

function retrieve(element::String; year::String="latest", type::String="scd", metastable::Bool=false, adas_type::String="adf11")
     println("Retrieving $adas_type - $type data for $element with year=$year and metastable:$metastable")
     database = retrieve_database(adas_type)
     data = retrieve_element(element,database)
     if !(data === nothing)
        return retrieve_data(data; year= year, type= type, metastable = metastable, adas_type = adas_type)
     else
        return nothing
     end
end

function retrieve_element(element, database)

            if !(lowercase(element) in collect(keys(database)))
                println("Cannot find element '$element' in database")
                println("Elements available are: ", collect(keys(database)))
                return nothing
            else
                return database[lowercase(element)]
            end
end

colors = Dict{}()
colors[0] = :blue
colors[1] = :red
colors[2] = :green
colors[3] = :magenta
function print_dict_keys(dict;level=0)
    for (k,v) in dict
        printstyled("   "^level,"- $k\n";color = colors[level])
        if typeof(v) <: Dict{Any,Any}
            print_dict_keys(v,level=level+1)
        end
    end
end

function show_database(element; adas_type="adf11")
    database = retrieve_database(adas_type)
    if !(lowercase(element) in collect(keys(database)))
        println("Cannot find element '$element' in database")
        println("Elements available are: ", collect(keys(database)))
        return
    else
        println("------ Element: $element --------")
        print_dict_keys(database[lowercase(element)])
        println("---------------------------------")
    end
end
function save_database(database,database_file)
    println("Saving database into $database_file")
    bson(database_file, database)
end

function load_database(path)
    println("Loading database from $path")
    return BSON.load(path)
end

function build_database(;kwargs...)
    build_adf11_database(kwargs...)
end

export retrieve_database
end



