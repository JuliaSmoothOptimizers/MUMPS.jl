using Libdl, LinearAlgebra, SparseArrays

using MPI

function __init__()
  _DEPS_FILE = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
  if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
  else
    error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
  end
  if (libmumps_simple == "julia_registryci_automerge")
    @error "MUMPS library not properly installed but module is loaded for AutoMerge"
  else
    global LIB_S = dlopen(libsmumps)
    global LIB_D = dlopen(libdmumps)
    global LIB_C = dlopen(libcmumps)
    global LIB_Z = dlopen(libzmumps)
  end
end

include("mumps_types.jl")
include("mumps_struc.jl")
include("interface.jl")
include("convenience.jl")
include("icntl_alibis.jl")
include("printing.jl")

include("exported_methods.jl")
