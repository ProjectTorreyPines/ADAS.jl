using Formatting, Printf, ADAS
"""
    Write a UEDGE ADAS file with the given Zmax and output stream.
    
    Parameters:
    - Zmax: Maximum atomic number (default is 100).
    - io: Output stream (default is stdout).
"""
function write_uedge_adas_file(filename::String="ADAS", args...; kw...)
    open(filename, "w") do io
        write_uedge_adas_file(io, args...; kw...)
    end
end
write_uedge_adas_file(filename::String, imp::Symbol; kw...) = write_uedge_adas_file(filename, ADAS.get_ionization_rate(imp; kw...), ADAS.get_recombination_rate(imp; kw...), ADAS.get_radiation_rate(imp; kw...), ADAS.get_cx_rate(imp; kw...))


function write_uedge_adas_file(filename::String, args...; kw...)
    open(filename, "w") do io
        write_uedge_adas_file(io, args...; kw...)
    end
end

function write_uedge_adas_file(io::Union{Base.TTY, IOStream},scd, acd, plt, ccd)
    Zmax = maximum(scd.Z) + 1
    @assert length(scd.Te) == length(acd.Te) == length(plt.Te) == length(ccd.Te)
    @assert length(scd.ne) == length(acd.ne) == length(plt.ne) == length(ccd.ne)
    Te = scd.Te
    ne = scd.ne
    n_ne = length(scd.ne)
    n_Te = length(scd.Te)
    scd_rates, acd_rates, plt_rates, cx_rates = make_rates(scd, acd, plt, ccd, Zmax, scd.Te, scd.ne)
    N = n_Te * n_ne * (Zmax+1)
    write_first_line(io)
    write_label(io, "$(basename(scd.scd.filepath)):$(basename(acd.acd.filepath)):$(basename(plt.plt.filepath)):$(basename(ccd.ccd.filepath))")
    write_header_dims(io)
    write_dims(io, n_ne, n_Te, Zmax)
    write_line_header(io, "real", Zmax + 1, "rtza")
    write_line_array(io, collect(0:Zmax))
    write_line_header(io, "real", Zmax + 1, "rtzn")
    write_line_array(io, collect(0:Zmax) .* 0 .+ Zmax)
    write_line_header(io, "real", Zmax + 1, "rtza2")
    write_line_array(io, collect(0:Zmax) .^ 2)

    write_line_header(io, "real", n_Te, "rtt")
    write_line_array(io, Te)
    write_line_header(io, "real", n_ne, "rtn")
    write_line_array(io, ne)
    write_line_header(io, "real", n_Te, "rtlt")
    write_line_array(io, Te; islog=true)
    write_line_header(io, "real", n_ne, "rtln")
    write_line_array(io, ne; islog=true)
    write_line_header(io, "real", N, "rtlsa")
    write_line_array(io, scd_rates; islog=true)
    write_line_header(io, "real", N, "rtlra")
    write_line_array(io, acd_rates; islog=true)
    write_line_header(io, "real", N, "rtlqa")
    write_line_array(io, plt_rates; islog=true)
    write_line_header(io, "real", N, "rtlcx")
    write_line_array(io, cx_rates; islog=true)
end


function write_line_header(io, b::AbstractString, n::Integer, c::AbstractString)
    # 2a8: two 8-character strings, i12: 12-digit integer, 4x: 4 spaces, a32: 32-character string
    @printf(io, "%-8.8s%-8.8s%12d    %-32.32s\n", "*cf:", b, n, c)
end

function write_line_array(io, arr; islog=false, line_max=5)
    # 2a8: two 8-character strings, i12: 12-digit integer, 4x: 4 spaces, a32: 32-character string
        c = 1
        for i in eachindex(arr)
        arr_ = arr[i]
        if islog
            arr_ = log(arr[i])
        end
        @printf(io," %18.13E", arr_)
        if c > line_max
            @printf(io, "\n")
            c = 1
        else
            c += 1
        end
    end
    if c> 1
        @printf(io, "\n")
    end
end

function make_rates(scd, acd, plt, ccd, Zmax, Te, ne)
    scd_rates = []
    for Z_ in 0:Zmax-1
        for ne_ in ne
            for Te_ in Te
                push!(scd_rates, scd(Z_, ne_, Te_))
            end
        end
    end

    for ne_ in ne
        for Te_ in Te
            push!(scd_rates, 1e-99)
        end
    end

    acd_rates = []
    for ne_ in ne
        for Te_ in Te
            push!(acd_rates, 1e-99)
        end
    end
    for Z_ in 1:Zmax
        for ne_ in ne
            for Te_ in Te
                push!(acd_rates, acd(Z_, ne_, Te_))
            end
        end
    end

    cx_rates = []
    for ne_ in ne
        for Te_ in Te
            push!(cx_rates, 1e-99)
        end
    end
    for Z_ in 1:Zmax
        for ne_ in ne
            for Te_ in Te
                push!(cx_rates, ccd(Z_, ne_, Te_) )
            end
        end
    end


    plt_rates = []
    for Z_ in 0:Zmax
        for ne_ in ne
            for Te_ in Te
                push!(plt_rates, plt(Z_, ne_, Te_) / ee)
            end
        end
    end
    return scd_rates, acd_rates, plt_rates, cx_rates
end



write_first_line(io) = print(io, "*cf:    char             120    label  \n")
write_label(io, label="") = print(io, "*cf:    make by J.Guterl | $label \n")
write_header_dims(io) = print(io,"*cf:    int                3    rtnt,rtnn,rtns\n ")
write_dims(io, n_ne, n_Te, Zmax) = print(io, "$(n_Te-1) $(n_ne-1) $(Zmax+1)\n")

x = 0.1
imp = :C
year = "latest"
scd = ADAS.get_ionization_rate(imp; year)
acd = ADAS.get_recombination_rate(imp; year)
plt = ADAS.get_radiation_rate(imp; year)
ccd = ADAS.get_cx_rate(imp; year)


io = stdout
scd_rates, acd_rates, plt_rates, cx_rates = make_rates(scd, acd, plt, ccd, Zmax, Te, ne)
@assert length(scd_rates) == length(acd_rates) == length(plt_rates) == length(cx_rates)
N = length(scd_rates)






using DelimitedFiles
a = DelimitedFiles.readdlm("examples/test.txt")
scd_test = DelimitedFiles.readdlm("examples/scd_test.txt")
te = DelimitedFiles.readdlm("examples/te_test.txt")
te_arr = Float64[]
for j in 1:7
for i in 1:6
    if te[j,i] isa Float64
        push!(te_arr, te[j,i])
    end
end
end
arr = Float64[]
for j in 1:814
for i in 1:6
    if a[j,i] isa Float64
        push!(arr, a[j,i])
    end
end
end

arr_scd = Float64[]
for j in 1:814
    for i in 1:6
        if scd_test[j, i] isa Float64
            push!(arr_scd, scd_test[j, i])
        end
    end
end


arr_ = reshape(arr, 41, 17, 7)
arr_scd = reshape(arr_scd, 41, 17, 7)
using Plots


plot()
for Z_ in 1:7
    plot!(te_arr, exp.(arr_scd[:, 7, Z_]); label="Z=$Z_", linestyle=:dash, linewidth=2.0, color=Z_)
    plot!(te_arr, scd.(Z_ .- 1 .+ 0 .* te_arr, 1e19 .+ 0 .* te_arr, te_arr); label="Z=$Z_", color=Z_)
end
plot!(; xscale=:log10)

plot()
for Z_ in 1:7
plot!(te_arr, exp.(arr_[:,7,Z_]), label="Z=$Z_", linestyle=:dash, linewidth=2.0, color=Z_)
plot!(te_arr, plt.(Z_ .-1 .+ 0 .* te_arr, 1e20 .+ 0 .* te_arr, te_arr) ./ (ee); label="Z=$Z_", color=Z_)
end
plot!(xscale=:log10)

using Plots

ee = 1.602176634e-19 # elementary charge in Coulombs
for Z_ in 1:6
end
plot!()