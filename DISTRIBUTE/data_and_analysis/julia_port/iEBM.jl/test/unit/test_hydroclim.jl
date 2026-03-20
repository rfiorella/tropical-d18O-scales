using Test
using iEBM

@testset "Hydroclim" begin
    @testset "compute_slope flat surface" begin
        lat = collect(Float64, -5:1:5)
        lon = collect(Float64, 0:1:10)
        elevation = zeros(length(lat), length(lon))

        slope = compute_slope(elevation, lat, lon)
        @test all(slope .== 0.0)
    end

    @testset "compute_slope uniform gradient" begin
        lat = collect(Float64, 0:1:10)
        lon = collect(Float64, 0:1:10)
        nlat = length(lat)
        nlon = length(lon)

        # Elevation increases linearly with latitude: elev = lat * 100
        elevation = zeros(nlat, nlon)
        for j in 1:nlon, i in 1:nlat
            elevation[i, j] = lat[i] * 100.0
        end

        slope = compute_slope(elevation, lat, lon)

        # Meridional gradient should be ~100 m/deg, zonal should be ~0
        # Slope = sqrt(100^2 + 0^2) = 100 at interior points
        @test slope[5, 5] ≈ 100.0 atol=1e-10
    end

    @testset "build_orog_mask" begin
        config = RunConfig()  # defaults: slope_threshold=400, elev_threshold=300
        nlat, nlon = 10, 10
        slope = fill(500.0, nlat, nlon)     # above threshold
        elevation = fill(400.0, nlat, nlon) # above threshold

        mask = build_orog_mask(slope, elevation, config)
        @test all(mask)

        # Below thresholds
        slope_low = fill(100.0, nlat, nlon)
        mask_low = build_orog_mask(slope_low, elevation, config)
        @test !any(mask_low)
    end

    @testset "partition_field" begin
        nlat, nlon = 5, 10
        field = ones(nlat, nlon)
        mask = falses(nlat, nlon)
        mask[3, 5] = true
        landfrac = ones(nlat, nlon)

        fixed, movable = partition_field(field, mask, landfrac)

        # Zonal mean is 1.0 everywhere, so anomaly is 0, fixed is 0
        @test all(fixed .== 0.0)
        @test all(movable .≈ field)
    end
end
