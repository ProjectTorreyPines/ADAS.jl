#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

using DataFrames
using Interpolations
using DelimitedFiles

function read_adf15(path::String; order::Int64=1)

    log10pec_dict = Dict{String, Dict{String, Any}}()

    # Open the file and read lines
    lines = readlines(path)
    header = split(lines[1])

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

        # Extract Z if possible
        Z = try
            parse(Int, join(header[4:end-3]))
        catch
            0
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
        for item in header[2:end]
            parts = split(item, '=')
            if length(parts) == 2
                header_dict[strip(parts[1])] = strip(parts[2])
            end
        end


        line_idx += 1

        # Densities
        dens = Float64[]
        while length(dens) < num_den
            dens = vcat(dens, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end
        dens = collect(dens)

        # Temperatures
        temp = Float64[]
        while length(temp) < num_temp
            temp = vcat(temp, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end
        temp = collect(temp)

        # PEC values
        PEC = Float64[]
        while length(PEC) < num_den * num_temp
            PEC = vcat(PEC, parse.(Float64, split(lines[line_idx])))
            line_idx += 1
        end
        PEC = reshape(PEC, num_den, num_temp) # reshape vector to matrix

        # Create PEC interpolator using log10 values
        # Define the spline type based on the value of `order` 
        # order = 1 -> Linear
        # order = 2 -> Quadratic
        # higher orders could be added (default order == 1) e.g. spline_type = order == 1 ? Linear() : order == 2 ? Quadratic() : Cubic()

        interp_type = order == 1 ? Linear() : Quadratic()

        # Interpolate log10 PEC with user-specified spline order
        pec_fun = interpolate(
            (log10.(dens), log10.(temp)),
            log10.(PEC),
            Gridded(interp_type),
        )

        # Populate dictionary
        log10pec_dict[string(header_dict["lam"])] = Dict(
            "log10 PEC fun" => pec_fun,
            "dens pnts" => dens,
            "temp pnts" => temp,
            "PEC pnts" => PEC,
            #"lambda [A]" => header_dict["lam"],
            "type" => lowercase(header_dict["type"]),
            "INDM" => parse(Int, get(header_dict, "INDM", "1"))
        )
    end

    return log10pe

end

function get_interpolated_value(log10pec_dict, lambda_input, dens_input, temp_input)
    # Step 1: Find the closest lambda in the dictionary keys
    lambda_keys = collect(parse.(Float64, keys(log10pec_dict)))

    closest_lambda = findmin(abs.(lambda_keys .- lambda_input))[2]
    lambda_key = string(lambda_keys[closest_lambda])

    # Step 2: Retrieve the dictionary entry for the closest lambda
    lambda_data = log10pec_dict[lambda_key]

    # Extract the interpolator and grid points from the entry
    interpolator = lambda_data["log10 PEC fun"]
    dens_points = lambda_data["dens pnts"]
    temp_points = lambda_data["temp pnts"]

    # Step 3: Interpolate
    # Ensure that dens_input and temp_input are within the range of grid points
    dens_input = clamp(dens_input, minimum(dens_points), maximum(dens_points))
    temp_input = clamp(temp_input, minimum(temp_points), maximum(temp_points))

    # Perform the interpolation
    interpolated_value = 10 .^ interpolator.(log10.(dens_input), log10.(temp_input))

    return interpolated_value
end


