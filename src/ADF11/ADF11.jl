#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#

include("types.jl")

include("parser.jl")

function vec2dict(adf11_files::Vector{<:adf11File})
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

# function get_adf11_database(;database_file = default_adf11_database_file, adf11_directory = default_adf11_directory)
#     if !isfile(database_file)
#         build_adf11_database(database_file = database_file, directory = adf11_directory)
#     end
#     return load_database(database_file)
# end

function retrieve_adf11_element_data(data::Dict; year::Union{String,Missing}="latest", type::String="scd", metastable::Bool=false)
    type = lowercase(type)
    types = collect(keys(adas_type_dict[:adf11]))
    if !(type in types)
        error("Cannot find type $type for adf11 format. Available types: $types")
    end
    if !(lowercase(type) in collect(keys(data)))
        error("Cannot find type '$(lowercase(type))' in available types: $(collect(keys(data)))")
    end

    data = data[type]
    if (year isa Missing) || (year == "latest")
        year = sort!(collect(keys(data)))[end]
    end

    if !(lowercase(year) in collect(keys(data)))
        println("Cannot find year '$year' in data")
        println("Year available are: ", collect(keys(data)))
        return nothing
    end

    data = data[year]

    if metastable
        if !("metastable" in collect(keys(data)))
            println("Cannot find metastable in data")
            return nothing
        else
            return data["metastable"]
        end
    else
        if !("ground" in collect(keys(data)))
            println("Cannot find ground in data")
            return nothing
        else
            return data["ground"]
        end
    end

end

get_data_type(data::adf11File{T}) where {T} = T
