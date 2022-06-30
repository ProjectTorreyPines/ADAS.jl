const adf11_types = Dict("acd" => "effective recombination",
"scd" => "effective ionization",
"prb" => "continuum radiation",
"plt" => "line radiation",
"ccd" => "thermal charge exchange",
"prc"=> "thermal charge exchange continuum radiation",
"pls"=> "line radiation in the SXR range",
"prs"=> "continuum radiation in the SXR range",
"brs"=> "continuum spectral bremstrahlung",
"fis"=> "sensitivity in the SXR range",
"pbs"=> "impurity bremsstrahlung in SXR range, also included in prs files")

function show_adf11_types()
    println("----------- adf11 types ------------")
    for (k,v) in adf11_types
        println("'$k' : $v")
    end
    print("--------------------------------------")
end

function read_adf11_file(filepath::String)
    
    if verbose
       println("Reading ADAS file: $filepath") 
    end

    @assert isfile(filepath) "File not found: $filepath"
    fid = open(filepath,"r") 
    lines = readlines(fid)
    close(fid)
    return lines
end


@redef struct ADASBlock
    header :: String
    content :: Vector{String}
end

@redef struct ADASHeader
    
    n_ions :: Int64
    n_ne :: Int64
    n_Te :: Int64
    imin :: Int64
    imax :: Int64
    details :: String
    function ADASHeader(block::ADASBlock)
        header = block.content[1]
        data = parse.(Int64,split(header)[1:5])
        details =  join(split(header)[6:end]," ")
        if verbose
            println("ADASHeader : data=$data; details=$details")
        end
        new(data...,details)
    end
end


function get_block_attr(block_header::String, attr::String; outputformat = Int64)
    rg = Regex("$attr\\s*=\\s*([\\d]+)")
    m = match(rg,block_header)
    if m === nothing
        return nothing
    else
        return parse.(outputformat,m.captures[1])
    end
end

function split_blocks(lines::Vector{String})
    blocks = Vector{ADASBlock}()
    iprevious = 1
    for i in 1:length(lines)
        if startswith(strip(lines[i]),"--") || startswith(strip(lines[i]),"C-") 
            push!(blocks, ADASBlock(lines[max(1,iprevious-1)],lines[iprevious:i-1]))
            iprevious = i+1
        end

    end
    if verbose; println("# of blocks = ", length(blocks)); end
    return blocks
end

function get_comments(lines::Vector{String})
    comments = Vector{String}()
    for i in 1:length(lines)
        if startswith(lines[i],"C") 
            push!(comments, lines[i])
        end
    end
    return comments
end

function get_units(comments::Vector{String})
    for c in comments
        m = match(Regex("IN UNITS OF"),c) # todo: capture directly units through regex
         if m != nothing
            return split(comments[3],"IN UNITS OF")[end]
        end
    end
    return ""
end


@redef struct ADASAxis
    log_Te :: Vector{Float64}
    log_ne :: Vector{Float64}
    Te :: Vector{Float64}
    ne :: Vector{Float64}
    unit_Te :: String
    unit_ne :: String
    function ADASAxis(block::ADASBlock, header :: ADASHeader)
    log_Te = Vector{Float64}()
    log_ne = Vector{Float64}()
    

    idx = 1
    while length(log_ne) < header.n_ne
        push!(log_ne,parse.(Float64,[n for n in split(block.content[idx])])...)
        idx += 1
    end
    while length(log_Te) < header.n_Te
        push!(log_Te,parse.(Float64,[t for t in split(block.content[idx])])...)
        idx += 1
    end

    new(log_Te, log_ne, 10 .^log_Te, 10 .^log_ne, "eV", "cm^-3")
    end
end

@redef struct ADASRates{T} 
    igrnd :: Union{Int64,Nothing}
    iptr :: Union{Int64,Nothing}
    Z :: Int64
    log_values :: T
    values :: T
    function ADASRates(block::ADASBlock, header::ADASHeader)
        Z = get_block_attr(block.header,"Z1")
        igrd = get_block_attr(block.header,"IGRD")
        iptr = get_block_attr(block.header,"IPRT")
        tmp = Vector{Float64}() #zeros(Float64,header.n_ne * header.n_Te)
        idx = 1
        while length(tmp) < header.n_ne * header.n_Te
            push!(tmp,parse.(Float64,[replace(r,"D" =>"e") for r in split(block.content[idx])])...)
            idx += 1
        end
        log_values = reshape(tmp,header.n_ne, header.n_Te)
        if verbose
             println("igrd = $igrd; iptr = $iptr; Z = $Z")
        end
        new{typeof(log_values)}(igrd,iptr,Z,log_values, 10.0 .^log_values)
    end
end

@redef struct adf11Data{T}
    filepath :: String
    header :: ADASHeader
    axis :: ADASAxis
    rates :: T
    units :: String
    function adf11Data(filepath::String; metastable=false)
        lines = read_adf11_file(filepath)
        comments = get_comments(lines)
        units = get_units(comments)
        blocks = split_blocks(lines)
        header = ADASHeader(blocks[1])
        if metastable
            axis = ADASAxis(blocks[3], header )
       

        rates = Dict{Int64,Vector{ADASRates}}()
   
        for block in blocks
            Z = get_block_attr(block.header,"Z1")
            if Z != nothing
                if !(Z in collect(keys(rates)))
                    rates[Z] = Vector{ADASRates}()
                end 
                push!(rates[Z],ADASRates(block,header)) 
            end
        end
    
        else
        axis = ADASAxis(blocks[2], header )
        rates = Dict{Int64,ADASRates}()
            
        for block in blocks    
            
            Z = get_block_attr(block.header,"Z1")
            if Z != nothing
                rates[Z] = ADASRates(block,header)
            end
        end
    end
        new{typeof(rates)}(filepath,header, axis, rates, units)
    end
end

@redef struct adf11File
    name :: String
    element :: String
    path :: String
    year :: String
    metastable :: Bool
    type ::  String
    data :: adf11Data
    function adf11File(filename::String, filepath::String)
        name = filename 
        element = match(Regex("(?<=\\_)(.*?)(?=\\.|\\#)"),filename).captures[1]
        type = match(Regex("(.*?)(?=[0-9]{2})"),filename).captures[1]
        year = match(Regex("([0-9]{2})"),filename).captures[1]
        if match(Regex("(?<=[0-9]{2})(r)(?=\\_)"),filename) !== nothing
            metastable = true
        else
            metastable = false
        end
        data = adf11Data(filepath; metastable = metastable)
        new(name, element,filepath,year,metastable,type,data)
    end
end


function get_list_adf11file(directory)
    adf11_files = Dict{String,adf11File}()
    if verbose; println("Parsing the directory $directory"); end
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if occursin(".dat",file)
                adf11_files[file] = adf11File(file, joinpath(root, file))
            end
    end
end
    return adf11_files
end


function generate_adf11_database(adf11_files :: Dict{String,adf11File})
database = Dict()
    for (name,file) in adf11_files
        if !(file.element in collect(keys(database)))
            database[file.element] = Dict()
        end
        if !(file.type in collect(keys(database[file.element])))
            database[file.element][file.type] = Dict()
        end
        if !(file.year in collect(keys(database[file.element][file.type])))
            database[file.element][file.type][file.year] = Dict()
        end
        if file.metastable
            database[file.element][file.type][file.year]["metastable"] = file
        else
            database[file.element][file.type][file.year]["ground"] = file
        end

    end  
    return database
end

function make_adf11_database(directory)
    print("Constructing database from directory: $directory")
    adf11_files = get_list_adf11file(directory)
    return generate_adf11_database(adf11_files)
end

function build_adf11_database(;database_file = default_adf11_database_file, directory = default_adf11_directory)
        database = make_adf11_database(directory)
        save_database(database, database_file)
end

function get_adf11_database(;database_file = default_adf11_database_file, adf11_directory = default_adf11_directory)
    if !isfile(database_file)
        build_adf11_database(database_file = database_file, directory = adf11_directory)
    end
    return load_database(database_file)
end

function retrieve_adf11_element_data(data; year::String="latest", type::String="scd", metastable::Bool=false)
    type = lowercase(type)
    types = collect(keys(adf11_types))
    if !(type in types)
        error("Cannot find type $type for adf11 format. Available types: $types")
    end
    if !(lowercase(type) in collect(keys(data)))
        println("Cannot find type '$element' in data")
        println("Types available are: ", collect(keys(data)))
        return nothing
    end

    data = data[type]
    if year == "latest"
        year = sort(collect(keys(data)))[end]
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