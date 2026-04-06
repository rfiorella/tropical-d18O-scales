using Test, MAT
using DeltapModel

@testset "Reference comparison" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures")

    for case_file in readdir(fixtures_dir; join=true)
        endswith(case_file, ".mat") || continue
        startswith(basename(case_file), "reference_") || continue

        @testset "$(basename(case_file))" begin
            data = matread(case_file)

            alpha_eq = data["alpha_eq"]

            deltap_bar, tau_bar = get_deltap_cpu(
                data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
                data["LAT"], data["LON"], data["LAT2"], data["LON2"],
                data["delta_e"], alpha_eq, data["Plim"], data["dmax"]
            )

            ref_dp  = data["deltap_bar"]
            ref_tau = data["tau_bar"]

            # Where ref is NaN, Julia should also be NaN OR have a valid value
            # (Julia may compute a few more valid points than Octave due to
            # minor interpolation differences at the P/Plim threshold)
            # But where Julia is NaN, ref should also be NaN:
            @test all(isnan.(ref_dp[isnan.(deltap_bar)]))
            @test all(isnan.(ref_tau[isnan.(tau_bar)]))

            # Compare where both are non-NaN
            valid_dp  = .!isnan.(ref_dp) .& .!isnan.(deltap_bar)
            valid_tau = .!isnan.(ref_tau) .& .!isnan.(tau_bar)

            n_valid_dp  = count(valid_dp)
            n_valid_tau = count(valid_tau)

            @test n_valid_dp > 0
            @test n_valid_tau > 0

            if n_valid_dp > 0
                max_err_dp = maximum(abs.(deltap_bar[valid_dp] .- ref_dp[valid_dp]))
                @test max_err_dp < 1e-4
                println("  deltap_bar max error: $max_err_dp ($n_valid_dp points)")
            end

            if n_valid_tau > 0
                max_err_tau = maximum(abs.(tau_bar[valid_tau] .- ref_tau[valid_tau]))
                @test max_err_tau < 1e-4
                println("  tau_bar max error: $max_err_tau ($n_valid_tau points)")
            end
        end
    end
end
