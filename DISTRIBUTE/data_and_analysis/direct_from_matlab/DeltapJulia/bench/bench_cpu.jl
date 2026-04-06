"""
Benchmark: naive CPU vs threaded/optimized CPU vs GPU kernel

Usage:
    julia --project bench/bench_cpu.jl
    julia --project -t4 bench/bench_cpu.jl    # 4 threads
"""

using MAT
using DeltapModel
using KernelAbstractions

fixtures_dir = joinpath(@__DIR__, "..", "test", "fixtures")

be, _ = select_backend()
println("Julia threads: ", Threads.nthreads())
println("GPU backend:   ", be)
println()

for fname in ["reference_tiny.mat", "reference_small.mat"]
    path = joinpath(fixtures_dir, fname)
    isfile(path) || continue

    data = matread(path)
    alpha_eq = data["alpha_eq"]
    args = (data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
            data["LAT"], data["LON"], data["LAT2"], data["LON2"],
            data["delta_e"], alpha_eq, data["Plim"], data["dmax"])

    A, B = size(data["LAT2"])
    println("=== $fname ($(A)×$(B) grid) ===")

    # Warmup
    get_deltap_cpu(args...)
    get_deltap_cpu_threaded(args...)
    run_gpu(args...; T=Float32, backend=CPU())

    # Benchmark naive CPU
    t_naive = @elapsed for _ in 1:3
        get_deltap_cpu(args...)
    end
    t_naive /= 3

    # Benchmark threaded CPU
    t_threaded = @elapsed for _ in 1:3
        get_deltap_cpu_threaded(args...)
    end
    t_threaded /= 3

    # Benchmark GPU kernel on CPU backend
    t_gpu_cpu = @elapsed for _ in 1:3
        run_gpu(args...; T=Float32, backend=CPU())
    end
    t_gpu_cpu /= 3

    println("  Naive CPU:       $(round(t_naive * 1000; digits=1)) ms")
    println("  Threaded CPU:    $(round(t_threaded * 1000; digits=1)) ms")
    println("  GPU kernel(CPU): $(round(t_gpu_cpu * 1000; digits=1)) ms")
    println("  Speedup (threaded/naive): $(round(t_naive / t_threaded; digits=2))×")

    # Benchmark on actual GPU if available
    if !(be isa KernelAbstractions.CPU)
        run_gpu(args...; T=Float32)  # warmup
        t_gpu = @elapsed for _ in 1:3
            run_gpu(args...; T=Float32)
        end
        t_gpu /= 3
        println("  GPU (Metal):     $(round(t_gpu * 1000; digits=1)) ms")
        println("  Speedup (GPU/naive): $(round(t_naive / t_gpu; digits=2))×")
    end

    println()
end
