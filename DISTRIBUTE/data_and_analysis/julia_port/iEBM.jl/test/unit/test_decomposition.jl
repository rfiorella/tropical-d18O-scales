using Test
using iEBM

@testset "Decomposition" begin
    @testset "spatial decompose: identity when old == new" begin
        # When old and new fields are identical, local and regional effects should be ~0
        n = 100
        dist = collect(range(0.0, 5000.0, length=n))
        tau = collect(range(0.0, 4.0, length=n))
        E0 = ones(n)
        nanindex = trues(n)

        local_threshold = 1000.0

        re, le = spatial_decompose(local_threshold,
            dist, tau, E0, nanindex,
            dist, tau, E0, nanindex)

        # Both effects should be ~0 when old == new
        @test abs(re) < 1e-10
        @test abs(le) < 1e-10
    end

    @testset "_taubar_simple" begin
        # Simple check: uniform E0, linear tau
        tau = collect(range(0.0, 3.0, length=100))
        E0 = ones(100)
        tb = iEBM._taubar_simple(tau, E0)
        # Should be positive and less than max tau
        @test tb > 0.0
        @test tb < 3.0
    end

    @testset "_integrate_tau" begin
        # With constant P/Fmag = mu, tau should be linear
        n = 50
        P = fill(1e-5, n)
        Fmag = fill(10.0, n)
        dx = 14.0
        tau = iEBM._integrate_tau(P, Fmag, dx)

        @test tau[1] == 0.0
        # tau[end] should be approximately (n-1) * dx * 1000 * P/Fmag
        expected = (n - 1) * dx * 1000 * 1e-5 / 10.0
        @test abs(tau[end] - expected) / expected < 1e-10
    end
end
