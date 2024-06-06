#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 runtest.jl (c) 2024
=#
using ADAS, Test
function test_basic_functionalities()
    data = retrieve_ADAS_data("C"; year="latest", type="scd", metastable=false)
    # show  available ADAS data
    show_ADAS_data()
    # show  available ADAS data for C
    show_ADAS_data(:C)
    # show adf11 format available
    show_adf11_types()
    zeff = ADAS.get_Zeff(:Kr)
    # retrieve cooling rate for W
    cr = ADAS.get_cooling_rates(:W)
    # retrieve abundance fraction
    af = ADAS.get_abundance_fraction(:Mo)
    true
end

function test_database_building()
    ADAS.build_ADAS_database(; parsed_data_path="test_parsed_data")
    true
end

@testset "Test basic functionalities of ADAS.jl" begin
    @test test_basic_functionalities()
    @test test_database_building()
end