#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#

module ADAS

using Logging
using MD5
using AbstractTrees
using FileIO
abstract type ADASRate end
abstract type ADASFile{T} end
struct ADASType{T} end

include("ADF11/ADF11.jl")
include("ADF15/ADF15.jl")
include("data.jl")
include("utils.jl")
include("plot_recipes.jl")

export ADASdata
export retrieve_ADAS_data, show_ADAS_data, build_ADAS_database, show_adf11_types

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
