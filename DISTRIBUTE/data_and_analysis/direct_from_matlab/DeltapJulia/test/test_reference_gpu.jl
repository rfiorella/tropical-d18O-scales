using Test, MAT
using DeltapModel
using KernelAbstractions

@testset "Reference comparison (GPU kernel, CPU backend)" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures")

    for case_file in readdir(fixtures_dir; join=true)
        endswith(case_file, ".mat") || continue
        startswith(basename(case_file), "reference_") || continue

        @testset "$(basename(case_file))" begin
            data = matread(case_file)
            alpha_eq = data["alpha_eq"]

            @testset "Float64" begin
                deltap_bar, tau_bar = run_gpu(
                    data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
                    data["LAT"], data["LON"], data["LAT2"], data["LON2"],
                    data["delta_e"], alpha_eq, data["Plim"], data["dmax"];
                    T=Float64, backend=CPU()
                )

                ref_dp  = data["deltap_bar"]
                ref_tau = data["tau_bar"]

                @test all(isnan.(ref_dp[isnan.(deltap_bar)]))
                @test all(isnan.(ref_tau[isnan.(tau_bar)]))

                valid_dp  = .!isnan.(ref_dp) .& .!isnan.(deltap_bar)
                valid_tau = .!isnan.(ref_tau) .& .!isnan.(tau_bar)

                @test count(valid_dp) > 0
                @test count(valid_tau) > 0

                if count(valid_dp) > 0
                    max_err = maximum(abs.(deltap_bar[valid_dp] .- ref_dp[valid_dp]))
                    @test max_err < 1e-4
                    println("  GPU(F64) deltap_bar max error: $max_err")
                end
                if count(valid_tau) > 0
                    max_err = maximum(abs.(tau_bar[valid_tau] .- ref_tau[valid_tau]))
                    @test max_err < 1e-4
                    println("  GPU(F64) tau_bar max error: $max_err")
                end
            end

            @testset "Float32" begin
                deltap_bar, tau_bar = run_gpu(
                    data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
                    data["LAT"], data["LON"], data["LAT2"], data["LON2"],
                    data["delta_e"], alpha_eq, data["Plim"], data["dmax"];
                    T=Float32, backend=CPU()
                )

                ref_dp  = data["deltap_bar"]
                ref_tau = data["tau_bar"]

                valid_dp  = .!isnan.(ref_dp) .& .!isnan.(deltap_bar)
                valid_tau = .!isnan.(ref_tau) .& .!isnan.(tau_bar)

                @test count(valid_dp) > 0
                @test count(valid_tau) > 0

                if count(valid_dp) > 0
                    max_err = maximum(abs.(deltap_bar[valid_dp] .- ref_dp[valid_dp]))
                    @test max_err < 1e-3  # Relaxed tolerance for Float32
                    println("  GPU(F32) deltap_bar max error: $max_err")
                end
                if count(valid_tau) > 0
                    max_err = maximum(abs.(tau_bar[valid_tau] .- ref_tau[valid_tau]))
                    @test max_err < 1e-3
                    println("  GPU(F32) tau_bar max error: $max_err")
                end
            end
        end
    end
end
