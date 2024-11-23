#=
Author: Luca Cappelli (cappellil@fusion.gat.com)
ADAS.jl (c) 2024
=#

using DataFrames
using Interpolations
using DelimitedFiles

path = "/home/cappellil/adf15_python/pec40#w_cl#w1.dat"
order = Int(1)

log10pec_dict = Dict{Int, Dict{String, Any}}()

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

    isel = parse(Int, header_dict["isel"])

    # Populate dictionary
    log10pec_dict[isel] = Dict(
        "log10 PEC fun" => pec_fun,
        "dens pnts" => dens,
        "temp pnts" => temp,
        "PEC pnts" => PEC,
        "lambda [A]" => header_dict["lam"],
        "type" => lowercase(header_dict["type"]),
        "INDM" => parse(Int, get(header_dict, "INDM", "1"))
    )
end

# Convert dictionary to DataFrame
out_data = [
    (
        "isel" => k,
        "lambda [A]" => v["lambda [A]"],
        "type" => v["type"],
        "log10 PEC fun" => v["log10 PEC fun"],
        "indm" => v["INDM"],
        "dens pnts" => v["dens pnts"],
        "temp pnts" => v["temp pnts"],
        "PEC pnts" => v["PEC pnts"]
    ) for (k, v) in log10pec_dict
]