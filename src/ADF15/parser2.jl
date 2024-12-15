#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

# version 1.0 :
# - coding with a few structures and datatypes, more julia-like coding should be implemented
# - only excitation is considered
# - no attributes are retrieved, nor transition 
# - data pahts defined in 'data.jl' inside variable 'data_directory'

function get_adf15_datafile(Element::Symbol, Z::Float64; bundling_model::String = "ic")

    element_string = lowercase(string(Element))

    path_data = data_directory[:adf15]
    
    pattern = ".*" * bundling_model * "#" * "[A-Za-z][a-z]?" * string(Int64(Z)) * ".dat"

    rg = Regex(pattern)

    # List all files in the directory
    files = readdir(path_data)

    # Filter files that match the regex
    matching_files = filter(f -> occursin(rg, f), files)

    # Print matching files
    if !isempty(matching_files)
        println("Matching files:")
        println(matching_files)

        fullname = path_data * "/" * matching_files[1]

    else
        println("No matching files found.")

        url = "https://open.adas.ac.uk/download/adf15/pec40]" *
        "[" * element_string * "/pec40]" *
        "[" * element_string * "_" * bundling_model * "]" *
        "[" * element_string * string(Int64(Z)) * ".dat"


        filename = "pec40#" * element_string * "_" * bundling_model * "#" * element_string * string(Int64(Z)) * ".dat"

        fullname = path_data * "/" * filename

        Downloads.download(url, fullname)

        println("Downloading: " * fullname)

    end

    return fullname
end

function read_adf15(path::String; order::Int64=1)

    pec_dict = Dict{String, Dict{String, Any}}()
    meta = Dict{String, Any}()

    # Open the file and read lines
    lines = readlines(path)
    header = split(lines[1])

    # Extract Z from header if possible
    Z = try
        parse(Int, join(header[4:end-3]))
    catch
        0
    end

    # Get the expected number of lines by reading the header
    num_lines = parse(Int, header[1])
    spec = strip(header[2], '/')

    line_idx = 2  # Start after the first header line

    for i in 1:num_lines # for each transition index (isel) skip lines you don't need to reach headers
        # Process the line with "isel" that is the header
        # Skip until "isel" appears in the line
        while !occursin("isel", lowercase(lines[line_idx]))
            line_idx += 1
        end

        while !occursin("isel", lowercase(lines[line_idx]))
            line_idx += 1
        end

        # header = ['wavelength (Angstrom) num_dens num_temp' , '...']
            
        # create header0 splitting header[1]

        # header[1] = 'wavelength (Angstrom) num_dens num_temp' 
            
        # header0 = [wavelength (Angstrom), num_dens, num_temp]


        header_dict = Dict{String, Any}()
        header = split(lines[line_idx], '/')
        header0 = split(header[1])

        header_dict["lam"] = parse(Float64, strip(header0[1], 'A'))
        num_den = parse(Int, header0[end-1])
        num_temp = parse(Int, header0[end])

        # Parse other parameters from header to header_dict
        # for instance emission type: excitation, Recombination, Charge Exchange
        for item in header[2:end]
            parts = split(item, '=')
            if length(parts) == 2
                header_dict[strip(parts[1])] = strip(parts[2])
            end
        end


        line_idx += 1

        # Densities (cm^-3)
        dens = Float64[]
        while length(dens) < num_den
            dens = vcat(dens, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end
        dens = collect(dens)

        # Temperatures (eV)
        temp = Float64[]
        while length(temp) < num_temp
            temp = vcat(temp, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end
        temp = collect(temp)

        # PEC values (cm^3 s^-1)
        PEC = Float64[]
        while length(PEC) < num_den * num_temp
            PEC = vcat(PEC, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end

        PEC = reshape(PEC, num_temp, num_den)' # reshape vector to matrix

        # Create PEC interpolator
        # Define the spline type based on the value of `order` 
        # order = 1 -> Linear
        # order = 2 -> Quadratic
        # higher orders could be added (default order == 1) e.g. spline_type = order == 1 ? Linear() : order == 2 ? Quadratic() : Cubic()

        interp_type = order == 1 ? Linear() : Quadratic()

        # Interpolate  PEC with user-specified spline order (to be developed)
        pec_fun = interpolate(
            (dens, temp),
            PEC,
            Gridded(interp_type),
        )
        
        pec_fun = interpolate((dens, temp), PEC, Gridded(interp_type))

        # Populate dictionary
        pec_dict[string(header_dict["lam"])] = Dict(
            "PEC fun" => pec_fun,
            "dens pnts" => dens,
            "temp pnts" => temp,
            "PEC pnts" => PEC,
            "type" => lowercase(header_dict["type"]),
            "INDM" => parse(Int, get(header_dict, "INDM", "1"))
        )

        meta["spec"] =  spec
        meta["Z"] = Z
        meta["num_den"] = num_den
        meta["num_temp"] = num_temp   

    end

    return pec_dict, meta
end


struct PhotonEmissivityCoefficients{U}
    pec::U
    Te::Vector{Float64}
    ne::Vector{Float64}
    wavelengths::Vector{Float64}
    Z::Float64
    imp::Symbol
end

function get_photon_emissivity_coeff(datafile::String; kw...)
    
    pec_dict, meta = read_adf15(datafile)
    wavelengths = parse.(Float64, keys(pec_dict))
    Nw = length(wavelengths)

    Z = meta["Z"]
    imp = Symbol(meta["spec"])

    ne = pec_dict[string(wavelengths[1])]["dens pnts"]
    Te = pec_dict[string(wavelengths[1])]["temp pnts"]

    itp_type = typeof(pec_dict[string(wavelengths[1])]["PEC fun"])

    pec = Vector{itp_type}(undef, Nw)

    for i = 1:Nw
        pec[i] = pec_dict[string(wavelengths[i])]["PEC fun"]
    end

    return PhotonEmissivityCoefficients{Vector{itp_type}}(pec, Te, ne, wavelengths, Z, imp)
end


function get_pec_interpolated_value(pec_dict, lambda_input, dens_input, temp_input)
    # Step 0: convert input to ADAS units
    dens_input = dens_input ./ 1e6 # conversion: m⁻³ -> cm⁻³ 
    
    # Step 1: Find the closest lambda in the dictionary keys
    lambda_keys = collect(parse.(Float64, keys(pec_dict)))

    closest_lambda = findmin(abs.(lambda_keys .- lambda_input))[2]
    lambda_key = string(lambda_keys[closest_lambda])

    # Step 2: Retrieve the dictionary entry for the closest lambda
    lambda_data = pec_dict[lambda_key]

    # Extract the interpolator and grid points from the entry
    interpolator = lambda_data["PEC fun"]
    dens_points = lambda_data["dens pnts"]
    temp_points = lambda_data["temp pnts"]

    # Step 3: Interpolate
    # Ensure that dens_input and temp_input are within the range of grid points
    dens_input = clamp_with_warning!(dens_input, dens_points, "dens_input")
    temp_input = clamp_with_warning!(temp_input, temp_points, "temp_input")

    # Perform the interpolation
    interpolated_value = interpolator.(dens_input, temp_input) # [cm³ s⁻¹]

    println("selected wavelength: " * "$lambda_key [Å]")

    return interpolated_value, lambda_key
end

function clamp_with_warning!(input, coordinate, variable_name)
    clamped_input = clamp.(input, minimum(coordinate), maximum(coordinate))
    if any(clamped_input .!= input)
        @warn "Values in $variable_name were outside the boundaries and have been clamped."
    end
    return clamped_input
end
