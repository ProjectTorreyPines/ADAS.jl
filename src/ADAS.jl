module ADAS


using Logging
using MD5
using AbstractTrees
using FileIO
abstract type ADASRate end
abstract type ADASFile{T} end
struct ADASType{T} end

include("ADF11/ADF11.jl")
include("data.jl")
include("cooling_rates.jl")
include("plot_recipes.jl")

export ADASdata
export retrieve_ADAS_data, show_ADAS_data, build_ADAS_database

end
