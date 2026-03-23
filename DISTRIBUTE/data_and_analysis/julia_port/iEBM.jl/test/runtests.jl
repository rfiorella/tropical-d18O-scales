using Test

@testset "iEBM Tests" begin
    @testset "Unit Tests" begin
        include("unit/test_config.jl")
        include("unit/test_interpolation.jl")
        include("unit/test_streamline_march.jl")
        include("unit/test_taubar_integration.jl")
        include("unit/test_terrestrial_efrac.jl")
        include("unit/test_decomposition.jl")
        include("unit/test_hydroclim.jl")
    end

    # Integration tests: run only when fixtures are available or explicitly requested
    if get(ENV, "IEBM_INTEGRATION", "0") == "1" ||
       isdir(joinpath(@__DIR__, "fixtures")) && !isempty(readdir(joinpath(@__DIR__, "fixtures")))
        @testset "Integration Tests" begin
            include("integration/test_houston_bbox.jl")
            include("integration/test_coord_list.jl")
            include("integration/test_decomposition_run.jl")
        end
    end
end
