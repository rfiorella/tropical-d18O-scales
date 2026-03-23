using Test
using iEBM

const FIXTURE_DIR = joinpath(@__DIR__, "..", "fixtures")
const COORD_INPUT = joinpath(FIXTURE_DIR, "houston_input.nc")
const COORD_FILE = joinpath(FIXTURE_DIR, "test_coords.csv")

@testset "Integration: Coordinate List" begin
    if !isfile(COORD_INPUT) || !isfile(COORD_FILE)
        @warn "Skipping coord-list integration test: fixtures not found"
        @test_skip true
        return
    end

    config = make_houston_config()
    # Override to use coord list mode
    # config = RunConfig(config; io=IOConfig(config.io; coord_file=COORD_FILE))

    @testset "runloop in coord_list mode" begin
        output, df_sl = runloop(config)
        @test output !== nothing
        @test df_sl !== nothing
    end
end
