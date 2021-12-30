using Libdl, LinearAlgebra, SparseArrays

using MPI


function __init__()
  if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
  else
    error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
  end
  global LIB_S = dlopen(libsmumps)
  global LIB_D = dlopen(libdmumps)
  global LIB_C = dlopen(libcmumps)
  global LIB_Z = dlopen(libzmumps)
end

include("mumps_types.jl")
include("mumps_struc.jl")
include("interface.jl")
include("convenience.jl")
include("icntl_alibis.jl")
include("printing.jl")

include("exported_methods.jl")
