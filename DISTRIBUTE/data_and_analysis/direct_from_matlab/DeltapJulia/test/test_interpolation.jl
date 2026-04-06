using Test
using DeltapModel

@testset "Bilinear interpolation" begin
    # 3×4 grid with known analytic values: f(lat, lon) = 2*lat + 3*lon
    lat_grid = [10.0, 20.0, 30.0]
    lon_grid = [100.0, 110.0, 120.0, 130.0]
    field = [2*lat + 3*lon for lat in lat_grid, lon in lon_grid]

    @testset "exact grid points" begin
        @test bilinear_interp(field, lat_grid, lon_grid, 10.0, 100.0) ≈ 320.0
        @test bilinear_interp(field, lat_grid, lon_grid, 20.0, 120.0) ≈ 400.0
        @test bilinear_interp(field, lat_grid, lon_grid, 30.0, 130.0) ≈ 450.0
    end

    @testset "midpoints" begin
        # f(15, 115) = 2*15 + 3*115 = 375
        @test bilinear_interp(field, lat_grid, lon_grid, 15.0, 115.0) ≈ 375.0
        # f(25, 125) = 2*25 + 3*125 = 425
        @test bilinear_interp(field, lat_grid, lon_grid, 25.0, 125.0) ≈ 425.0
    end

    @testset "arbitrary interior point" begin
        # f(12, 108) = 2*12 + 3*108 = 348
        @test bilinear_interp(field, lat_grid, lon_grid, 12.0, 108.0) ≈ 348.0
    end

    @testset "edge clamping" begin
        # Query below grid minimum — should clamp
        v = bilinear_interp(field, lat_grid, lon_grid, 5.0, 95.0)
        @test v ≈ field[1, 1]  # clamped to corner
    end
end

@testset "RegularGrid bilinear interpolation" begin
    lat_grid = [10.0, 20.0, 30.0]
    lon_grid = [100.0, 110.0, 120.0, 130.0]
    field = [2*lat + 3*lon for lat in lat_grid, lon in lon_grid]
    grid = RegularGrid(lat_grid, lon_grid)

    @test bilinear_interp(field, grid, 10.0, 100.0) ≈ 320.0
    @test bilinear_interp(field, grid, 20.0, 120.0) ≈ 400.0
    @test bilinear_interp(field, grid, 15.0, 115.0) ≈ 375.0
    @test bilinear_interp(field, grid, 25.0, 125.0) ≈ 425.0
    @test bilinear_interp(field, grid, 12.0, 108.0) ≈ 348.0
    @test bilinear_interp(field, grid, 5.0, 95.0) ≈ field[1, 1]  # clamped
end

@testset "RegularLookup 1D interpolation" begin
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 2.0, 4.0, 6.0]
    lut = RegularLookup(x, y)

    @test DeltapModel.interp1_extrap(lut, 0.5) ≈ 1.0
    @test DeltapModel.interp1_extrap(lut, 1.5) ≈ 3.0
    @test DeltapModel.interp1_extrap(lut, -1.0) ≈ -2.0
    @test DeltapModel.interp1_extrap(lut, 4.0) ≈ 8.0
end

@testset "1D interpolation with extrapolation" begin
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 2.0, 4.0, 6.0]  # y = 2x

    @test DeltapModel.interp1_extrap(x, y, 0.5) ≈ 1.0
    @test DeltapModel.interp1_extrap(x, y, 1.5) ≈ 3.0
    @test DeltapModel.interp1_extrap(x, y, -1.0) ≈ -2.0  # extrapolate below
    @test DeltapModel.interp1_extrap(x, y, 4.0) ≈ 8.0    # extrapolate above

    # Vectorized version
    result = DeltapModel.interp1_extrap(x, y, [0.0, 1.5, 3.0])
    @test result ≈ [0.0, 3.0, 6.0]
end
