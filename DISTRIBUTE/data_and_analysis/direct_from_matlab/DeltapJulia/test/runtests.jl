using Test

@testset "DeltapModel" begin
    include("test_interpolation.jl")

    # Reference tests require .mat fixtures from Octave
    fixtures_dir = joinpath(@__DIR__, "fixtures")
    ref_files = filter(f -> startswith(f, "reference_") && endswith(f, ".mat"),
                       readdir(fixtures_dir))
    if !isempty(ref_files)
        include("test_reference.jl")
        include("test_reference_threaded.jl")
        include("test_reference_gpu.jl")
    else
        @warn "No reference .mat files found in test/fixtures/ — skipping reference tests"
    end
end
