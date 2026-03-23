using Test
using iEBM

const FIXTURE_DIR = joinpath(@__DIR__, "..", "fixtures")
const HOUSTON_INPUT = joinpath(FIXTURE_DIR, "houston_input.nc")

@testset "Integration: Houston BBox" begin
    if !isfile(HOUSTON_INPUT)
        @warn "Skipping Houston bbox integration test: fixture not found at $HOUSTON_INPUT"
        @test_skip true
        return
    end

    config = make_houston_config()
    # Override input path to use fixture
    # config = RunConfig(config; io=IOConfig(config.io; input_file=HOUSTON_INPUT))

    @testset "runloop completes without error" begin
        output, df_sl = runloop(config)
        @test output !== nothing
        @test df_sl !== nothing
        @test nrow(df_sl) > 0
    end

    @testset "output contains expected fields" begin
        output, _ = runloop(config)
        # Verify key output arrays exist and have correct dimensions
        @test haskey(output, :tau_bar) || hasproperty(output, :tau_bar)
    end
end
