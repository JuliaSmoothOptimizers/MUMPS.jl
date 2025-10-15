# Test utilities for MUMPS.jl
# Provides helpers to silence verbose C/MUMPS output during tests.

using MUMPS

"""quiet_icntl() -> Vector{Int}
Return an ICNTL vector configured to suppress most MUMPS printing.
This uses `get_icntl(verbose=false)` and forces icntl[1:4]=0, icntl[4]=0.
"""
function quiet_icntl()
  icntl = get_icntl(verbose = false)
  # ensure output streams are suppressed
  icntl[1:4] .= 0
  return icntl
end

"""quiet_mumps(T; sym=mumps_unsymmetric)
Create a `Mumps{T}` instance using quiet ICNTL settings for tests.
"""
function quiet_mumps(::Type{T}; sym = mumps_unsymmetric) where {T}
  icntl = quiet_icntl()
  cntl = T <: Float32 ? default_cntl32 : default_cntl64
  return Mumps{T}(sym, icntl, cntl)
end
