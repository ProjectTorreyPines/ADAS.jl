module ADAS

using BSON
using Logging
import Term.Trees: Tree
abstract type ADASFile{T} end
struct ADASType{T} end 
export ADASdata
export retrieve_ADAS_data

include("ADF11/ADF11.jl")
include("data.jl")

end



