using Test
using iEBM

@testset "Interpolation" begin
    @testset "bilinear on analytic function" begin
        # Create a grid and evaluate sin(lat_rad) * cos(lon_rad)
        lat = collect(Float64, -90:1:90)
        lon = collect(Float64, 0:1:359)
        nlat = length(lat)
        nlon = length(lon)

        data = zeros(Float64, nlat, nlon)
        for j in 1:nlon, i in 1:nlat
            data[i, j] = sind(lat[i]) * cosd(lon[j])
        end

        interp = RegularGridInterp(lat, lon, data)

        # Test at grid points (should be exact)
        @test interp(0.0, 0.0) ≈ sind(0.0) * cosd(0.0) atol=1e-12
        @test interp(45.0, 90.0) ≈ sind(45.0) * cosd(90.0) atol=1e-12
        @test interp(-30.0, 180.0) ≈ sind(-30.0) * cosd(180.0) atol=1e-12

        # Test at midpoints (bilinear will have some error on nonlinear functions,
        # but should be close on a 1-degree grid)
        val = interp(30.5, 45.5)
        exact = sind(30.5) * cosd(45.5)
        @test abs(val - exact) < 1e-4
    end

    @testset "longitude wraparound" begin
        lat = collect(Float64, -90:1:90)
        lon = collect(Float64, 0:1:359)
        nlat = length(lat)
        nlon = length(lon)

        # Constant field = 1.0
        data = ones(Float64, nlat, nlon)
        interp = RegularGridInterp(lat, lon, data)

        # Should interpolate to 1.0 everywhere, including near boundaries
        @test interp(0.0, 0.5) ≈ 1.0 atol=1e-12
        @test interp(0.0, 359.5) ≈ 1.0 atol=1e-12
    end

    @testset "build_interpolators" begin
        lat = collect(Float64, -10:1:10)
        lon = collect(Float64, 0:1:10)
        n = length(lat)
        m = length(lon)

        E = ones(Float64, n, m)
        P = fill(2.0, n, m)
        UQ = fill(3.0, n, m)
        VQ = fill(4.0, n, m)
        lfrac = fill(0.5, n, m)

        Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
            lat, lon, E, P, UQ, VQ, lfrac)

        @test Efit(0.0, 5.0) ≈ 1.0 atol=1e-12
        @test Pfit(0.0, 5.0) ≈ 2.0 atol=1e-12
        @test Fmagfit(0.0, 5.0) ≈ 5.0 atol=1e-12  # sqrt(9+16)
    end
end
