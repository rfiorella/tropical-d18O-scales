using Test
using iEBM

@testset "Streamline March" begin
    @testset "uniform zonal flow" begin
        # Set up a uniform zonal flow field at the equator
        # UQ = 10 kg/m/s, VQ = 0 -> pure eastward transport
        # P = 1e-5 kg/m^2/s (uniform), E = 1e-5 kg/m^2/s
        # Fmag = 10 kg/m/s
        lat = collect(Float64, -10:0.5:10)
        lon = collect(Float64, 0:0.5:359.5)
        nlat = length(lat)
        nlon = length(lon)

        P_val = 1e-5
        E_val = 1e-5
        UQ_val = 10.0
        VQ_val = 0.001  # small to avoid division by zero
        Fmag_val = sqrt(UQ_val^2 + VQ_val^2)

        P = fill(P_val, nlat, nlon)
        E = fill(E_val, nlat, nlon)
        UQ = fill(UQ_val, nlat, nlon)
        VQ = fill(VQ_val, nlat, nlon)
        lfrac = ones(nlat, nlon)  # all land

        Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
            lat, lon, E, P, UQ, VQ, lfrac)

        dx = 14.0   # km
        Dx = dx / 111.0
        taumax = 8.0
        Nmax = ceil(Int, 25000.0 / dx)

        result = StreamlineResult(Nmax)
        march_streamline!(result, 0.0, 180.0, Nmax, taumax, Dx, dx,
                          Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

        # tau should increase monotonically
        nsteps = result.nsteps
        @test nsteps > 10
        @test all(diff(result.tau[1:nsteps]) .>= 0)

        # Check that tau integration is approximately:
        # dtau/ds ≈ P/Fmag, so tau ≈ (P/Fmag) * dist_m
        # For first few steps, tau should be close to analytic
        mu_analytic = P_val / Fmag_val
        # After n steps at dx km each, tau ≈ mu * n * dx * 1000
        n_check = 10
        tau_expected = mu_analytic * n_check * dx * 1000
        @test abs(result.tau[n_check + 1] - tau_expected) / tau_expected < 0.01
    end

    @testset "StreamlineResult reset" begin
        result = StreamlineResult(100)
        result.tau[1] = 1.0
        result.nsteps = 50
        reset!(result)
        @test result.tau[1] == 0.0
        @test result.nsteps == 0
    end
end
