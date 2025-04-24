function dump_data(element, directory, data)
    file_path = get_parsed_data_filepath(element, directory)
    println("dumping data for element $element in $file_path ")
    FileIO.save(file_path, data)
    return nothing
end

function dump_data(directory, data)
    for (el, d) in data
        dump_data(el, directory, d)
    end
    return nothing
end

get_parsed_data_filepath(element, directory) = joinpath(directory, "$element.jld2")

function load_data(path::String)
    @debug "Loading jld2 file: $path"
    return FileIO.load(path) #BSON.load(path, @__MODULE__)
end

load_data(element::String, path::String) = load_data(get_parsed_data_filepath(element, path))

function read_file(filepath::String)
    @debug "Reading file: $filepath"
    @assert isfile(filepath) "File not found: $filepath"
    fid = open(filepath, "r")
    lines = readlines(fid)
    close(fid)
    return lines
end



function get_file_list(directory, FileType::Type)::Vector{FileType}
    println("walking and looking for file of type $FileType in directory: $directory")
    files_ = Vector{FileType}()
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if occursin(".dat", file)
                push!(files_, FileType(joinpath(root, file)))
            end
        end
    end
    return files_
end

function _get_file_list(directory)
    files_ = Vector{String}()
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if occursin(".dat", file)
                push!(files_, file)
            end
        end
    end
    return files_
end

function get_adas_filetype(type::Symbol)
    type == :adf15 && return adf15File
    type == :adf11 && return adf11File
    error("Cannot find ADASFileType for type $type")
end

get_file_list(adas_type::Symbol) = get_file_list(data_paths[:raw_data][adas_type], get_adas_filetype(adas_type))

get_file_list(directory::String, adas_type::Symbol) = get_file_list(directory, get_adas_filetype(adas_type))

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
