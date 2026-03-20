using Test
using iEBM

@testset "Config" begin
    @testset "default construction" begin
        config = RunConfig()
        @test config.io.run_name == "houston"
        @test config.grid.deg_per_lat == 0.25
        @test config.isotope.streamline_dx_km == 14.0
        @test config.isotope.streamline_max_tau == 8.0
        @test config.decomp.local_threshold_km == 1000.0
        @test config.constants.Lv == 2.54e6
        @test config.constants.g == 9.81
    end

    @testset "make_houston_config" begin
        config = make_houston_config()
        @test config.isotope.compute_method == "bbox"
        @test config.isotope.bbox_lat_range == [28.0, 32.0]
        @test config.isotope.bbox_lon_range == [263.0, 267.0]
    end

    @testset "validation passes for defaults" begin
        config = RunConfig()
        @test validate(config) === nothing
    end

    @testset "validation catches bad compute_method" begin
        config = RunConfig(isotope=IsotopeConfig(compute_method="invalid"))
        @test_throws ErrorException validate(config)
    end

    @testset "validation catches bad dtau_fraction" begin
        config = RunConfig(decomp=DecompConfig(dtau_fraction=1.5))
        @test_throws ErrorException validate(config)
    end

    @testset "keyword construction" begin
        config = RunConfig(
            io = IOConfig(run_name="test_run"),
            isotope = IsotopeConfig(streamline_dx_km=10.0)
        )
        @test config.io.run_name == "test_run"
        @test config.isotope.streamline_dx_km == 10.0
        @test config.grid.deg_per_lat == 0.25  # unchanged default
    end
end
