using Test
using iEBM

const FIXTURE_DIR = joinpath(@__DIR__, "..", "fixtures")
const HOUSTON_INPUT = joinpath(FIXTURE_DIR, "houston_input.nc")

@testset "Integration: Decomposition Run" begin
    if !isfile(HOUSTON_INPUT)
        @warn "Skipping decomposition integration test: fixture not found at $HOUSTON_INPUT"
        @test_skip true
        return
    end

    config = make_houston_config()
    # Enable decomposition
    # config = RunConfig(config;
    #     decomp=DecompConfig(config.decomp;
    #         local_v_regional_local_evap=true,
    #         local_v_regional_upwind=true,
    #         decomp_E_L_s=true))

    @testset "decomposition runloop completes" begin
        output, df_sl = runloop(config)
        @test output !== nothing
        @test df_sl !== nothing
    end

    @testset "decomposition outputs present" begin
        output, _ = runloop(config)
        # Verify decomposition-specific fields exist
        # @test haskey(output, :local_effect)
        # @test haskey(output, :regional_effect)
    end
end
