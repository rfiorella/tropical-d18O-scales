using Test
using iEBM

@testset "Terrestrial E Fraction" begin
    @testset "all-land streamline" begin
        n = 100
        # Uniform exponentially decaying wp
        tau = collect(range(0.0, 4.0, length=n))
        E0 = ones(n)
        nanindex = trues(n)

        # wp = E * exp(-tau) / trapz(E * exp(-tau)), matching Python normalization
        wp_raw = E0 .* exp.(-tau)
        wp_sum = iEBM._trapz(wp_raw)
        wp = wp_raw ./ wp_sum

        # coast_step = n-1 matches Python semantics: in Python, coast_step for an
        # all-land streamline of n steps (0-based 0..n-1) is set to step = n-1.
        # Julia's march_streamline! now also produces n-1 for this case.
        coast_step = n - 1
        land_check = 1.0

        tau_wtd, land_frac, E_terr = terrestrial_efrac(
            land_check, wp, coast_step, tau, E0, nanindex, 10)

        @test land_frac ≈ 1.0 atol=0.02
        # E_terr uses sum(wp[1:coast_step]) which can slightly differ from 1.0
        # because wp is trapz-normalized but E_terr uses plain sum (matches Python)
        @test E_terr > 0.9
        @test E_terr < 1.1
        @test tau_wtd > 0.0
    end

    @testset "coast_step via march_streamline" begin
        # End-to-end test: all-land uniform flow, verify coast_step semantics
        lat = collect(Float64, -10:0.5:10)
        lon = collect(Float64, 0:0.5:359.5)
        nlat = length(lat)
        nlon = length(lon)

        P = fill(1e-5, nlat, nlon)
        E = fill(1e-5, nlat, nlon)
        UQ = fill(10.0, nlat, nlon)
        VQ = fill(0.001, nlat, nlon)
        lfrac = ones(nlat, nlon)  # all land

        Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
            lat, lon, E, P, UQ, VQ, lfrac)

        dx = 14.0
        Dx = dx / 111.0
        Nmax = 200
        result = StreamlineResult(Nmax)
        march_streamline!(result, 0.0, 180.0, Nmax, 8.0, Dx, dx,
                          Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

        # All land -> coast_step should be set to nsteps - 1
        # (matching Python's 0-based semantics where coast_step = step after loop)
        @test result.coast_step == result.nsteps - 1
        @test result.coast_step > 0

        # wp[1:coast_step] should not include more elements than exist
        @test result.coast_step <= length(result.wp)
    end

    @testset "ocean point returns zeros" begin
        n = 50
        tau = collect(range(0.0, 2.0, length=n))
        E0 = ones(n)
        nanindex = trues(n)
        wp = fill(1.0 / n, n)

        tau_wtd, land_frac, E_terr = terrestrial_efrac(
            0.0, wp, 0, tau, E0, nanindex, 5)

        @test tau_wtd == 0.0
        @test land_frac == 0.0
    end
end
