using Test

@testset "iEBM Tests" begin
    include("unit/test_config.jl")
    include("unit/test_interpolation.jl")
    include("unit/test_streamline_march.jl")
    include("unit/test_taubar_integration.jl")
    include("unit/test_terrestrial_efrac.jl")
    include("unit/test_decomposition.jl")
    include("unit/test_hydroclim.jl")
end
