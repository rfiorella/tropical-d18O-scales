"""
Stage 4 profiling: benchmark all backends on 0.5° global grid (720×360).

Usage:
    julia --project bench/bench_halfdeg.jl
    julia --project -t8 bench/bench_halfdeg.jl    # 8 threads
"""

using MAT
using DeltapModel
using KernelAbstractions

fixtures_dir = joinpath(@__DIR__, "..", "test", "fixtures")
ref_file = joinpath(fixtures_dir, "reference_halfdeg.mat")

if !isfile(ref_file)
    @error "reference_halfdeg.mat not found — run scripts/generate_halfdeg_testcase.m first"
    exit(1)
end

data = matread(ref_file)
alpha_eq = data["alpha_eq"]
args = (data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
        data["LAT"], data["LON"], data["LAT2"], data["LON2"],
        data["delta_e"], alpha_eq, data["Plim"], 20000)

A, B = size(data["LAT2"])
be, _ = select_backend()

println("="^60)
println("Stage 4 Profiling: 0.5° global grid ($(A)×$(B) = $(A*B) points)")
println("Julia threads: ", Threads.nthreads())
println("GPU backend:   ", be)
println("dmax:          ", 20000, " km")
println("Nmax:          ", ceil(Int, 20000 / 15.0))
println("="^60)
println()

# --- Validate against Octave reference ---
println("--- Validation ---")
ref_dp = data["deltap_bar"]
ref_tau = data["tau_bar"]
n_ref_valid = count(.!isnan.(ref_dp))
println("Octave reference: $n_ref_valid / $(length(ref_dp)) valid points")
println()

# --- Benchmark functions ---
function bench(name, f; n_runs=1)
    print("  $name: warming up... ")
    f()  # warmup
    GC.gc()
    println("benchmarking ($n_runs runs)...")
    times = Float64[]
    for r in 1:n_runs
        t = @elapsed result = f()
        push!(times, t)
    end
    t_mean = sum(times) / length(times)
    t_min = minimum(times)
    println("    mean: $(round(t_mean; digits=2))s, min: $(round(t_min; digits=2))s")
    return t_min
end

function validate_result(name, deltap_bar, tau_bar, ref_dp, ref_tau)
    valid_dp = .!isnan.(ref_dp) .& .!isnan.(deltap_bar)
    n = count(valid_dp)
    if n > 0
        max_err = maximum(abs.(deltap_bar[valid_dp] .- ref_dp[valid_dp]))
        println("    $name deltap_bar: max_err=$max_err ($n valid points)")
        return max_err
    else
        println("    $name: no valid points to compare!")
        return NaN
    end
end

# --- Naive CPU ---
println("--- Naive CPU ---")
t_naive = bench("Naive CPU", () -> get_deltap_cpu(args...))
dp_naive, tau_naive = get_deltap_cpu(args...)
validate_result("Naive CPU", dp_naive, tau_naive, ref_dp, ref_tau)
println()

# --- Threaded CPU ---
println("--- Threaded CPU ($(Threads.nthreads()) threads) ---")
t_threaded = bench("Threaded CPU", () -> get_deltap_cpu_threaded(args...))
dp_thr, tau_thr = get_deltap_cpu_threaded(args...)
validate_result("Threaded CPU", dp_thr, tau_thr, ref_dp, ref_tau)
println()

# --- GPU kernel on CPU backend ---
println("--- GPU kernel (CPU backend, Float32) ---")
t_gpu_cpu = bench("GPU(CPU,F32)", () -> run_gpu(args...; T=Float32, backend=CPU()))
dp_gcpu, tau_gcpu = run_gpu(args...; T=Float32, backend=CPU())
validate_result("GPU(CPU,F32)", dp_gcpu, tau_gcpu, ref_dp, ref_tau)
println()

# --- GPU kernel on CPU backend, Float64 ---
println("--- GPU kernel (CPU backend, Float64) ---")
t_gpu_cpu64 = bench("GPU(CPU,F64)", () -> run_gpu(args...; T=Float64, backend=CPU()))
dp_gcpu64, tau_gcpu64 = run_gpu(args...; T=Float64, backend=CPU())
validate_result("GPU(CPU,F64)", dp_gcpu64, tau_gcpu64, ref_dp, ref_tau)
println()

# --- Metal GPU ---
if !(be isa KernelAbstractions.CPU)
    println("--- Metal GPU (Float32) ---")
    t_metal = bench("Metal(F32)", () -> run_gpu(args...; T=Float32); n_runs=3)
    dp_metal, tau_metal = run_gpu(args...; T=Float32)
    validate_result("Metal(F32)", dp_metal, tau_metal, ref_dp, ref_tau)
    println()

    # Try different workgroup sizes
    println("--- Metal GPU workgroup size tuning ---")
    for wg in [32, 64, 128, 256, 512]
        try
            run_gpu(args...; T=Float32, workgroup_size=wg)  # warmup
            t = @elapsed run_gpu(args...; T=Float32, workgroup_size=wg)
            println("  workgroup=$wg: $(round(t; digits=3))s")
        catch e
            println("  workgroup=$wg: FAILED ($e)")
        end
    end
    println()
end

# --- Fixed-step comparison (nsteps = 1500, no tau<10 early exit) ---
println("="^60)
println("FIXED-STEP COMPARISON (nsteps=1500, tau_stop=Inf)")
println("="^60)
println()

NSTEPS_FIXED = 1500

println("--- Naive CPU (fixed 1500 steps) ---")
t_naive_fx = bench("Naive CPU (fixed)", () -> get_deltap_cpu(args...; fixed_nsteps=NSTEPS_FIXED))
dp_naive_fx, _ = get_deltap_cpu(args...; fixed_nsteps=NSTEPS_FIXED)
validate_result("Naive CPU (fixed)", dp_naive_fx, dp_naive_fx, ref_dp, ref_tau)
println()

println("--- Threaded CPU ($(Threads.nthreads()) threads, fixed 1500 steps) ---")
t_threaded_fx = bench("Threaded CPU (fixed)", () -> get_deltap_cpu_threaded(args...; fixed_nsteps=NSTEPS_FIXED))
dp_thr_fx, _ = get_deltap_cpu_threaded(args...; fixed_nsteps=NSTEPS_FIXED)
validate_result("Threaded CPU (fixed)", dp_thr_fx, dp_thr_fx, ref_dp, ref_tau)
println()

println("--- GPU kernel (CPU backend, Float32, fixed 1500 steps) ---")
t_gpu_cpu_fx = bench("GPU(CPU,F32,fixed)", () -> run_gpu(args...; T=Float32, backend=CPU(), fixed_nsteps=NSTEPS_FIXED))
dp_gcpu_fx, _ = run_gpu(args...; T=Float32, backend=CPU(), fixed_nsteps=NSTEPS_FIXED)
validate_result("GPU(CPU,F32,fixed)", dp_gcpu_fx, dp_gcpu_fx, ref_dp, ref_tau)
println()

if !(be isa KernelAbstractions.CPU)
    println("--- Metal GPU (Float32, fixed 1500 steps) ---")
    t_metal_fx = bench("Metal(F32,fixed)", () -> run_gpu(args...; T=Float32, fixed_nsteps=NSTEPS_FIXED); n_runs=3)
    dp_metal_fx, _ = run_gpu(args...; T=Float32, fixed_nsteps=NSTEPS_FIXED)
    validate_result("Metal(F32,fixed)", dp_metal_fx, dp_metal_fx, ref_dp, ref_tau)
    println()
end

println("--- Fixed-step summary (vs adaptive) ---")
println("  Naive CPU:     $(round(t_naive; digits=2))s adaptive → $(round(t_naive_fx; digits=2))s fixed  ($(round(t_naive_fx/t_naive; digits=2))× ratio)")
println("  Threaded CPU:  $(round(t_threaded; digits=2))s adaptive → $(round(t_threaded_fx; digits=2))s fixed  ($(round(t_threaded_fx/t_threaded; digits=2))× ratio)")
println("  GPU(CPU,F32):  $(round(t_gpu_cpu; digits=2))s adaptive → $(round(t_gpu_cpu_fx; digits=2))s fixed  ($(round(t_gpu_cpu_fx/t_gpu_cpu; digits=2))× ratio)")
if !(be isa KernelAbstractions.CPU)
    println("  Metal GPU:     $(round(t_metal; digits=2))s adaptive → $(round(t_metal_fx; digits=2))s fixed  ($(round(t_metal_fx/t_metal; digits=2))× ratio)")
end
println()

# --- Summary ---
println("="^60)
println("SUMMARY")
println("="^60)
println("  Naive CPU:       $(round(t_naive; digits=2))s")
println("  Threaded CPU:    $(round(t_threaded; digits=2))s  ($(round(t_naive/t_threaded; digits=1))× speedup)")
println("  GPU kernel(CPU): $(round(t_gpu_cpu; digits=2))s  ($(round(t_naive/t_gpu_cpu; digits=1))× speedup)")
if !(be isa KernelAbstractions.CPU)
    println("  Metal GPU:       $(round(t_metal; digits=2))s  ($(round(t_naive/t_metal; digits=1))× speedup)")
end
println()
println("Grid points:  $(A*B)")
pts_per_sec_naive = A*B / t_naive
println("Throughput (naive):    $(round(Int, pts_per_sec_naive)) points/s")
if !(be isa KernelAbstractions.CPU)
    pts_per_sec_metal = A*B / t_metal
    println("Throughput (Metal):    $(round(Int, pts_per_sec_metal)) points/s")
end
