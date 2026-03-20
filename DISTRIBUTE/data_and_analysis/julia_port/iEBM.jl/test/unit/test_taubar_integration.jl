using Test
using iEBM

@testset "TauBar Integration" begin
    @testset "compute_taubar_at_point with uniform mu" begin
        # With uniform mu and E0, tau_bar should equal:
        # integral(tau * exp(-tau) dtau) / integral(exp(-tau) dtau)
        # which for tau in [0, T] is known analytically
        n = 1000
        tau_max = 5.0
        tausteps = collect(range(0.0, tau_max, length=n))
        mu_taugrid = ones(n)      # uniform mu = 1
        E0_taugrid = ones(n)      # uniform E0 = 1

        taubar = compute_taubar_at_point(mu_taugrid, E0_taugrid, tausteps)

        # Analytic: integral(tau * exp(-tau), 0, T) / integral(exp(-tau), 0, T)
        # = [1 - (1+T)*exp(-T)] / [1 - exp(-T)]
        T = tau_max
        analytic = (1 - (1 + T) * exp(-T)) / (1 - exp(-T))
        @test abs(taubar - analytic) / analytic < 1e-3
    end

    @testset "dist_to_tau_coords" begin
        n = 50
        tau = collect(range(0.0, 4.0, length=n))
        nanindex = trues(n)
        mu = fill(0.001, n)
        wp = fill(1.0 / n, n)
        E0 = ones(n)

        mu_tg, wp_tg, E0_tg, tausteps = dist_to_tau_coords(tau, nanindex, mu, wp, E0)

        @test length(mu_tg) == n
        @test length(tausteps) == n
        @test tausteps[1] ≈ 0.0 atol=1e-12
        @test tausteps[end] ≈ 4.0 atol=1e-12
    end

    @testset "linear_idx_to_degrees" begin
        deg_grid = collect(Float64, -90:1:90)
        # Index 0 -> -90, index 181 -> 90
        result = linear_idx_to_degrees([0.0, 90.5, 181.0], deg_grid)
        @test result[1] ≈ -90.0 atol=1e-10
        @test result[end] ≈ 90.0 atol=1e-10
    end
end
