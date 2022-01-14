using Libdl, LinearAlgebra, SparseArrays

using MPI

function __init__()
  _DEPS_FILE = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
  if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
  else
    haskey(ENV, "MUMPS_PREFIX") && error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
  end
end

if haskey(ENV, "MUMPS_PREFIX")
  const libsmumps = joinpath(ENV["MUMPS_PREFIX"], "lib/libsmumps.$dlext")
  const libdmumps = joinpath(ENV["MUMPS_PREFIX"], "lib/libdmumps.$dlext")
  const libcmumps = joinpath(ENV["MUMPS_PREFIX"], "lib/libcmumps.$dlext")
  const libzmumps = joinpath(ENV["MUMPS_PREFIX"], "lib/libzmumps.$dlext")
else
  using MUMPS_jll
end

include("mumps_types.jl")
include("mumps_struc.jl")
include("interface.jl")
include("convenience.jl")
include("icntl_alibis.jl")
include("printing.jl")

include("exported_methods.jl")
