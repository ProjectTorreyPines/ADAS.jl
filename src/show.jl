
Base.show(io::IO, data::ADASDatabase) = print(io, "ADAS database")
Base.show(io::IO, ::MIME"text/plain", db::ADASDatabase) = AbstractTrees.print_tree(db)


AbstractTrees.printnode(io::IO, file::adf11File) = printstyled(io, basename(file.filepath); bold=true)
AbstractTrees.printnode(io::IO, file::adf15File) = printstyled(io, basename(file.filepath); bold=true)
AbstractTrees.printnode(io::IO, ::Dict{String,adf11File}) = printstyled(io,"┐")
AbstractTrees.printnode(io::IO, ::Dict{String,Any}) = printstyled(io,"┐")

AbstractTrees.printnode(io::IO, ::Dict{String,Dict{String,adf11File}}) = printstyled(io,"┐")
AbstractTrees.printnode(io::IO, ::Dict{String,Dict{String,Dict{String,ADAS.adf11File}}}) = printstyled(io,"┐")

Base.show(io::IO, data::adf11File) = print(io, "ADAS data | adf11 : $(basename(data.filepath)) ")
Base.show(io::IO, data::adf15File) = print(io, "ADAS data | adf15 : $(basename(data.filepath)) ")
Base.show(io::IO, ::MIME"text/plain", data::adf11File) = print(io, "ADAS data | adf11 : $(basename(data.filepath)) ")
Base.show(io::IO, ::MIME"text/plain", data::adf15File) = print(io, "ADAS data | adf15 : $(basename(data.filepath)) ")

function show_ADAS_data(element; adas_type=:adf11, kw...) 
    AbstractTrees.print_tree(ADASdb(element, adas_type); maxdepth=5, kw...)
end
show_ADAS_files(; adas_type=:adf11, kw...) = AbstractTrees.print_tree(_get_file_list(data_paths[:raw_data][adas_type]); kw...)
