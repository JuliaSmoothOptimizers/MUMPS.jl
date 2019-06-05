export default_icntl, default_cntl32, default_cntl64, Mumps, get_icntl,
       finalize, factorize!, solve!, solve,
       associate_matrix!, associate_rhs!, get_solution,
       mumps_unsymmetric, mumps_definite, mumps_symmetric,
       MUMPSException

using LinearAlgebra
using MPI
using SparseArrays

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
end


mumps_lib = :libmumps_simple
macro mumps_call(func, args...)
  quote
    ccall(($func, $mumps_lib), $(args...))
  end |> esc
end


"Exception type raised in case of error."
mutable struct MUMPSException <: Exception
  msg :: AbstractString
end

const MUMPSValueDataType = Union{Float32, Float64, ComplexF32, ComplexF64}
const MUMPSIntDataType   = Int32


# See MUMPS User's Manual Section 5.1.
"Default integer parameters."
default_icntl = zeros(Int32, 40);
default_icntl[1]  =  6;  # Output stream for error messages
default_icntl[2]  =  0;  # Output stream for diagonstics/stats/warnings
default_icntl[3]  =  6;  # Output stream for global info on host
default_icntl[4]  =  2;  # Output level for errors/warnings/diagnostics
default_icntl[5]  =  0;  # 0 = assembled matrix, 1 = elemental format
default_icntl[6]  =  7;  # permutation/scaling in analysis (7 = automatic)
default_icntl[7]  =  7;  # pivot order for factorization (7 = automatic)
default_icntl[8]  = 77;  # scaling in analysis/factorization (77 = automatic)
default_icntl[9]  =  1;  # 1: solve Ax=b, otherwise A'x=b
default_icntl[10] =  0;  # max number of iterative refinement steps
default_icntl[11] =  0;  # > 0: return stats collected during solve
default_icntl[12] =  0;  # ordering during analysis
default_icntl[13] =  0;  # >0: do not use ScaLapack on root frontal matrix
default_icntl[14] = 20;  # % workspace increase during analysis/fact
default_icntl[15] =  0;  # (not used)
default_icntl[16] =  0;  # (not used)
default_icntl[17] =  0;  # (not used)
default_icntl[18] =  0;  # 0 = matrix assembled on host
default_icntl[19] =  0;  # 1 = return Schur complement on host
default_icntl[20] =  0;  # 0 = dense rhs, 1 = sparse rhs
default_icntl[21] =  0;  # 0 = solution overwrites rhs, 1 = keep distributed
default_icntl[22] =  0;  # 0 = in core, 1 = out of core
default_icntl[23] =  0;  # max working memory
default_icntl[24] =  0;  # 0: null pivot=error
default_icntl[25] =  0;  # -1: compute nullspace basis
default_icntl[26] =  0;  # condense rhs on Schur variables (see 19)
default_icntl[27] = -8;  # blocking size for multiple rhs (<0: value * (-2))
default_icntl[28] =  0;  # 1: sequential analysis, 2: parallel, 0: automatic
default_icntl[29] =  0;  # ordering for parallel analysis (see 28)
default_icntl[30] =  0;  # compute entries of the inverse
default_icntl[31] =  0;  # discard factors after factorization (can't solve!)
default_icntl[32] =  0;  # (not used)
default_icntl[33] =  0;  # compute determinant
default_icntl[34] =  0;  # (not used)
default_icntl[35] =  0;  # (not used)
default_icntl[36] =  0;  # (not used)
default_icntl[37] =  0;  # (not used)
default_icntl[38] =  0;  # (not used)
default_icntl[39] =  0;  # (not used)
default_icntl[40] =  0;  # (not used)

# See MUMPS User's Manual Section 5.2.
# icntl[1] will be set to its default value if left at -1.
"Default single precision real parameters"
default_cntl32 = zeros(Float32, 15);
default_cntl32[1] = -1;    # relative threshold for numerical pivoting
default_cntl32[2] = sqrt(eps(Float32));  # tolerance for iterative refinement
default_cntl32[3] =  0.0;  # threshold to detect null pivots
default_cntl32[4] = -1.0;  # threshold for static pivoting (<0: disable)
default_cntl32[5] =  0.0;  # what null pivots are reset to
# default_cntl32[6-15] are not used.

"Default double precision real parameters"
default_cntl64 = zeros(Float64, 15);
default_cntl64[1] = -1;    # relative threshold for numerical pivoting
default_cntl64[2] = sqrt(eps(Float64));  # tolerance for iterative refinement
default_cntl64[3] =  0.0;  # threshold to detect null pivots
default_cntl64[4] = -1.0;  # threshold for static pivoting (<0: disable)
default_cntl64[5] =  0.0;  # what null pivots are reset to
# default_cntl64[6-15] are not used.

# Symbols for symmetry
"""Constant indicating that a general unsymmetric matrix will be
analyzed and factorized"""
const mumps_unsymmetric = 0;

"""Constant indicating that a symmetric definite matrix will be
analyzed and factorized"""
const mumps_definite    = 1;

"""Constant indicating that a general symmetric matrix will be
analyzed and factorized"""
const mumps_symmetric   = 2;


"""Abstract type representing a factorization with MUMPS.
All constructor arguments are optional. By default a general
unsymmetric matrix will be analyzed/factorized with default
integer and real parameters"""
mutable struct Mumps{Tv <: MUMPSValueDataType}
  __id    :: Int               # Pointer to MUMPS struct as an Int. Do not touch.
  __sym   :: Int32             # Value of sym used by Mumps.
  icntl   :: Array{Int32,1}    # Integer control parameters.
  cntl    :: Union{Array{Float32,1}, Array{Float64,1}}    # Real control parameters.
  n       :: Int32             # Order of the matrix factorized.
  infog   :: Array{Int32,1}
  rinfog  :: Union{Array{Float32,1}, Array{Float64,1}}
  nnz     :: Int               # Number of nonzeros in factors.
  det     :: Tv
  err     :: Int

  function Mumps{Tv}(sym :: Int,
                     icntl :: Array{Int32,1},
                     cntl  :: Union{Array{Float32,1}, Array{Float64,1}}) where {Tv <: MUMPSValueDataType}

    MPI.Initialized() || throw(MUMPSException("Initialize MPI first"));

    # Set default pivot threshold if required.
    if cntl[1] == -1
      cntl[1] = (sym == mumps_definite) ? 0.0 : 0.01
    end

    if Tv == Float32
      id = @mumps_call(:mumps_initialize_float, Ptr{Nothing},
                       (Int32, Ptr{Int32}, Ptr{Float32}), sym, icntl, cntl);
      rinfog = zeros(Float32, 20);
    elseif Tv == Float64
      id = @mumps_call(:mumps_initialize_double, Ptr{Nothing},
                       (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);
      rinfog = zeros(Float64, 20);
    elseif Tv == ComplexF32
      id = @mumps_call(:mumps_initialize_float_complex, Ptr{Nothing},
                       (Int32, Ptr{Int32}, Ptr{Float32}), sym, icntl, cntl);
      rinfog = zeros(Float32, 20);
    else
      id = @mumps_call(:mumps_initialize_double_complex, Ptr{Nothing},
                       (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);
      rinfog = zeros(Float64, 20);
    end

    id == C_NULL && throw(MUMPSException("Error allocating MUMPS structure"));

    infog = zeros(Int32, 40);

    id = reinterpret(Int, id)
    self = new{Tv}(id, Int32(sym), icntl, cntl, 0, infog, rinfog, 0, Tv(0), 0);
    finalizer(finalize, self)  # Destructor.
    return self;
  end
end


"Obtain an array of integer control parameters."
function get_icntl(;
                   det :: Bool=false,       # Compute determinant.
                   verbose :: Bool=false,   # Output intermediate info.
                   ooc :: Bool=false,       # Store factors out of core.
                   itref :: Int=0,          # Max steps of iterative refinement.
                   )
  icntl = default_icntl[:];
  icntl[33] = det ? 1 : 0;
  if !verbose
    icntl[1:4] .= 0;
  end
  icntl[22] = ooc ? 1 : 0;
  icntl[10] = itref;
  return icntl;
end


import Base.finalize

for (fname, elty) in ((:mumps_finalize_float, Float32),
                      (:mumps_finalize_double, Float64),
                      (:mumps_finalize_float_complex, ComplexF32),
                      (:mumps_finalize_double_complex, ComplexF64))

  @eval begin

    "Terminate a Mumps instance."
    function finalize(mumps :: Mumps{$elty})
      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call($(string(fname)), Nothing, (Ptr{Nothing},), id)
      mumps.__id = C_NULL;
      return mumps
    end

  end
end


include("MUMPS_lib.jl")
