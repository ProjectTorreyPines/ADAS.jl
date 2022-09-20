

const parsed_data_directory = Dict(:adf11 => joinpath(@__DIR__,"../parsed_data/adf11"))
const data_directory = Dict(:adf11 => joinpath(@__DIR__,"../data/adf11"))
mutable struct ADASData
    adf11 :: ADF11
    paths :: Dict{String,Dict{Symbol,String}}
end

const ADASdata = ADASData(Dict(),Dict("parsed_data" => parsed_data_directory, "raw_data" => data_directory))

(db :: ADASData)(element; adas_type::Symbol  = :adf11) = (db :: ADASData)(element, adas_type)
(db :: ADASData)(element, adas_type::String) = (db :: ADASData)(element, Symbol(adas_type))

function (ADASdata :: ADASData)(element, adas_type::Symbol)
    @assert adas_type == :adf11 "getting data of type $adas_type is not implemented yet..."
    retrieve_element_data(element, getfield(ADASdata, adas_type), adas_type)
end

function retrieve_element_data(element::String, data , adas_type::Symbol)
    element = lowercase(element)
    @debug "looking for $element in $(keys(data))"
    
    if element ∉ keys(data)
        @debug "typeof(get_element_data(element,ADASdata.paths,adas_type)) = $(typeof(get_element_data(element,ADASdata.paths,adas_type)))"
        data[element] = get_element_data(element,ADASdata.paths,adas_type)
    end
    @assert element ∈ keys(data) "Cannot find element '$element' in database. Elements available are: $(collect(keys(data)))"
    return data[element]
end



function retrieve_element_data(args...; adas_type = :adf11, kwargs...)
    if adas_type == :adf11
        return retrieve_adf11_element_data(args...;kwargs...)
    else
        error("getting data of type $adas_type is not implemented yet...")
    end
end
retrieve_ADAS_data(element::String; year::String="latest", type::String="scd", metastable::Bool=false, adas_type::Symbol=:adf11) = retrieve_element_data(ADASdata(element,adas_type); year= year, type= type, metastable = metastable, adas_type = adas_type)

function get_element_data(element,paths,adas_type)
    filepath = get_data_filepath( element,paths["parsed_data"][adas_type])
    @debug "Getting element data from file: $filepath"
    if !isfile(filepath)
        make_database(element,paths,adas_type)
    end

    return load_data(filepath)
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

function show_ADAS_data(element; adas_type=:adf11)
    data = ADASdata(element,adas_type)
    t = Tree(data ,title="ADAS data: $element",
    title_style="magenta",
    guides_style="yellow")
    print(t)
end

function show_ADAS_data(; adas_type=:adf11)
    for element in keys(get_database(ADASdata.paths["raw_data"][adas_type],adas_type))
        show_ADAS_data(element; adas_type=adas_type)
    end
end



function dump_data(element, directory, data)
        file_path = get_data_filepath(element,directory)
        bson(file_path, data)
    end

function dump_data(directory, data)
    for (el,d) in data
        dump_data(el,directory,d)
    end
end

get_data_filepath(element, directory) = joinpath(directory,"$element.bson")

function load_data(path::String)
    @debug "Loading bson file $path"
    return BSON.load(path, @__MODULE__)
end

load_data(element, path::String) = load_data(get_data_filepath(element,path))

function read_adas_file(filepath::String)
    @debug "Reading ADAS file: $filepath"
    @assert isfile(filepath) "File not found: $filepath"
    fid = open(filepath,"r") 
    lines = readlines(fid)
    close(fid)
    return lines
end

function get_list_file(directory, FileType)
    files_ = Vector{FileType}()
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if occursin(".dat",file)
                push!(files_,FileType(file, joinpath(root, file)))
            end
    end
end
    return files_
end

function make_database(paths, adas_type) 
    data = get_database(paths["raw_data"][adas_type],adas_type)
    dump_data(paths["parsed_data"][adas_type], data)
end

function make_database(element, paths, adas_type) 
    element = lowercase(element)
    data = get_database(paths["raw_data"][adas_type], adas_type)
    @assert element ∈ keys(data)
    dump_data(element, paths["parsed_data"][adas_type], data)
end

function build_database()
    for adas_type in [f for f in fieldnames(ADASData) if f != :paths]
        build_database(ADASdata.paths,adas_type)
    end
end

function build_database(paths,adas_type)
        println("buidling database for $adas_type files")
        data = get_database(paths["raw_data"][adas_type],adas_type)
        dump_data(paths["parsed_data"][adas_type], data)
end
get_database(directory::String, adas_type) = get_database(directory::String, ADASType{adas_type}())
get_database(directory::String, ::ADASType{:adf11}) = get_database(get_list_file(directory,adf11File))

function Base.show(io::IO, ::MIME"text/plain", data::ADASData)
    print(io, "ADAS data: adf11: $(collect(keys(data.adf11)))")
end

function Base.show(io::IO, data::ADASData)
    print(io, "ADAS data: adf11: $(join(collect(keys(data.adf11)),"; "))")
end

function Base.show(io::IO, ::MIME"text/plain", file::adf11File)
    println(io," ADAS adf11 data: $(file.name)")
    for field in [f for f in fieldnames(adf11File) if f != :data]
    println(io," "^1*"└─ $(field): $(getfield(file,field))")
    end
end

function Base.show(io::IO, data::adf11File)
    print(io, "ADAS data | adf11 : $(data.name) ")
end