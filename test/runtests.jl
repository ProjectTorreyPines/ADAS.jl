#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
 Company: General Atomics
 runtest.jl (c) 2024
=#

using ADAS, Test

function test_basic_functionalities()
    data = retrieve_ADAS_data("C"; year="latest", type="scd", metastable=false)
    # show  available ADAS data
    show_ADAS_files()
    show_ADAS_files(; adas_type=:adf15)

    # show  available ADAS data for C
    show_ADAS_data(:C)

    # show adf11 format available
    show_adf11_types()
    zeff = ADAS.get_Zeff(:Kr)

    # retrieve cooling rate for W
    cr = ADAS.get_cooling_rates(:W)

    # retrieve abundance fraction
    af = ADAS.get_abundance_fraction(:Mo)

    # retrieve emission rate
    e = ADAS.get_emission_rate(:cr; Z=0, model="llu");
    return true
end

function test_database_building()
    set_adas_parsed_data_directory!("test_parsed_data")
    ADAS.build_ADAS_database()
    return true
end

@testset "Test basic functionalities of ADAS.jl" begin
    @test test_basic_functionalities()
    @test test_database_building()
end