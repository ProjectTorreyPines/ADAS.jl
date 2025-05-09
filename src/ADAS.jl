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
using Downloads
using Interpolations

abstract type ADASRate end
abstract type ADASFileType{T} end

abstract type ADASFile{T} end
struct ADASType{T} end
abstract type ADASRates{T} end
abstract type ADASHeader end 
abstract type ADASBlock end

const adas_type_dict = Dict()

include("ADF11/ADF11.jl")
include("ADF15/ADF15.jl")
include("io.jl")
include("data.jl")
include("utils.jl")
include("plot_recipes.jl")
include("show.jl")

export ADASdata
export retrieve_ADAS_data, show_ADAS_data, build_ADAS_database, show_adf11_types
export set_adas_parsed_data_directory!, set_adas_raw_data_directory!
const document = Dict() #export list of all names in module
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
