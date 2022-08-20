using Libdl, LinearAlgebra, SparseArrays

using MPI

if haskey(ENV, "JULIA_MUMPS_LIBRARY_PATH")
  println("Custom Installation")
  const libsmumps = joinpath(ENV["JULIA_MUMPS_LIBRARY_PATH"], "libsmumps.$dlext")
  const libdmumps = joinpath(ENV["JULIA_MUMPS_LIBRARY_PATH"], "libdmumps.$dlext")
  const libcmumps = joinpath(ENV["JULIA_MUMPS_LIBRARY_PATH"], "libcmumps.$dlext")
  const libzmumps = joinpath(ENV["JULIA_MUMPS_LIBRARY_PATH"], "libzmumps.$dlext")
  const libmumps_common = joinpath(ENV["JULIA_MUMPS_LIBRARY_PATH"], "libmumps_common.$dlext")
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
