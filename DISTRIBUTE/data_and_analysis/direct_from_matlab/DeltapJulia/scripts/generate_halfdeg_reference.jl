"""Generate Julia CPU reference for the 0.5° synthetic case."""

using MAT, DelimitedFiles

# Activate project
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DeltapModel

fixtures_dir = joinpath(@__DIR__, "..", "test", "fixtures")
input_file = joinpath(fixtures_dir, "synthetic_halfdeg_inputs.mat")
data = matread(input_file)

# Load alpha_eq from CSV
csv_path = joinpath(@__DIR__, "..", "data", "alpha_eq.csv")
alpha_eq_raw = readdlm(csv_path, ',', Float64; skipstart=1)
println("Input grid: ", size(data["LAT2"]))
println("alpha_eq: ", size(alpha_eq_raw))

println("Running naive CPU (this is the reference)...")
t = @elapsed deltap_bar, tau_bar = get_deltap_cpu(
    data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
    data["LAT"], data["LON"], data["LAT2"], data["LON2"],
    data["delta_e"], alpha_eq_raw, data["Plim"], 20000
)

valid = .!isnan.(deltap_bar)
println("Elapsed: $(round(t; digits=1))s")
println("Valid points: $(count(valid)) / $(length(deltap_bar))")
if count(valid) > 0
    println("deltap_bar range: [$(minimum(deltap_bar[valid])), $(maximum(deltap_bar[valid]))]")
end

# Save as reference
ref_file = joinpath(fixtures_dir, "reference_halfdeg.mat")
matwrite(ref_file, Dict(
    "E" => data["E"], "P" => data["P"], "UQ" => data["UQ"], "VQ" => data["VQ"],
    "Tcond" => data["Tcond"], "LAT" => data["LAT"], "LON" => data["LON"],
    "LAT2" => data["LAT2"], "LON2" => data["LON2"],
    "delta_e" => data["delta_e"], "alpha_eq" => alpha_eq_raw,
    "Plim" => data["Plim"], "dmax" => data["dmax"],
    "deltap_bar" => deltap_bar, "tau_bar" => tau_bar
))
println("Saved $ref_file")
