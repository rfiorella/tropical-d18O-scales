module DeltapModel

using DelimitedFiles
using Statistics
import Metal

include("interpolation.jl")
include("inpaint.jl")
include("types.jl")
include("streamline_cpu.jl")
include("streamline_threaded.jl")
include("streamline_gpu.jl")
include("gpu_driver.jl")

export bilinear_interp, inpaint_nans!, get_deltap_cpu, get_deltap_cpu_threaded,
       RegularGrid, RegularLookup, DeltapFields, run_gpu, select_backend

end
