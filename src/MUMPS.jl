module MUMPS

export default_icntl, default_cntl32, default_cntl64, Mumps, get_icntl,
       finalize, factorize, solve,
       associate_matrix, associate_rhs, associate_rhs!, get_solution,
       mumps_unsymmetric, mumps_definite, mumps_symmetric,
       MUMPSException

using MPI

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
end

# libjmumps.dylib should be on your LD_LIBRARY_PATH.
mumps_lib = "libmumps_simple";
macro mumps_call(func, args...)
  quote
    ccall(($func, $mumps_lib), $(args...))
  end
end


"Exception type raised in case of error."
type MUMPSException <: Exception
  msg :: AbstractString
end

typealias MUMPSValueDataType Union{Float32, Float64, Complex64, Complex128};
typealias MUMPSIntDataType   Union{Int64};


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
mumps_unsymmetric = 0;

"""Constant indicating that a symmetric definite matrix will be
analyzed and factorized"""
mumps_definite    = 1;

"""Constant indicating that a general symmetric matrix will be
analyzed and factorized"""
mumps_symmetric   = 2;


"""Abstract type representing a factorization with MUMPS.
All constructor arguments are optional. By default a general
unsymmetric matrix will be analyzed/factorized with default
integer and real parameters"""
type Mumps{Tv <: MUMPSValueDataType}
  __id    :: Int               # Pointer to MUMPS struct as an Int. Do not touch.
  __sym   :: Int32             # Value of sym used by Mumps.
  icntl   :: Array{Int32,1}    # Integer control parameters.
  cntl    :: Union{Array{Float32,1}, Array{Float64,1}}    # Real control parameters.
  n       :: Int32             # Order of the matrix factorized.
  infog   :: Array{Int32,1}
  rinfog  :: Union{Array{Float32,1}, Array{Float64,1}}
  nnz     :: Int               # Number of nonzeros in factors.
  det     :: AbstractFloat
  err     :: Int

  function Mumps(sym :: Int,
                 icntl :: Array{Int32,1},
                 cntl  :: Union{Array{Float32,1}, Array{Float64,1}})

    MPI.Initialized() || throw(MUMPSException("Initialize MPI first"));

    # Set default pivot threshold if required.
    if cntl[1] == -1
      cntl[1] = (sym == mumps_definite) ? 0.0 : 0.01
    end

    if Tv == Float32
      id = @mumps_call(:mumps_initialize_float, Ptr{Void},
                       (Int32, Ptr{Int32}, Ptr{Float32}), sym, icntl, cntl);
      rinfog = zeros(Float32, 20);
    elseif Tv == Float64
      id = @mumps_call(:mumps_initialize_double, Ptr{Void},
                       (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);
      rinfog = zeros(Float64, 20);
    elseif Tv == Complex64
      id = @mumps_call(:mumps_initialize_float_complex, Ptr{Void},
                       (Int32, Ptr{Int32}, Ptr{Float32}), sym, icntl, cntl);
      rinfog = zeros(Float32, 20);
    else
      id = @mumps_call(:mumps_initialize_double_complex, Ptr{Void},
                       (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);
      rinfog = zeros(Float64, 20);
    end

    id == C_NULL && throw(MUMPSException("Error allocating MUMPS structure"));

    infog = zeros(Int32, 40);

    id = reinterpret(Int, id)
    self = new(id, int32(sym), icntl, cntl, 0, infog, rinfog, 0, 0, 0);
    finalizer(self, finalize);  # Destructor.
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
    icntl[1:4] = 0;
  end
  icntl[22] = ooc ? 1 : 0;
  icntl[10] = itref;
  return icntl;
end


import Base.finalize

"Terminate a Mumps instance."
function finalize{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv})
  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    @mumps_call(:mumps_finalize_float, Void, (Ptr{Void},), id);
  elseif Tv == Float64
    @mumps_call(:mumps_finalize_double, Void, (Ptr{Void},), id);
  elseif Tv == Complex64
    @mumps_call(:mumps_finalize_float_complex, Void, (Ptr{Void},), id);
  else
    @mumps_call(:mumps_finalize_double_complex, Void, (Ptr{Void},), id);
  end
  mumps.__id = C_NULL;
  return mumps
end


"""Register the matrix `A` with the `Mumps` object `mumps`.
This function makes it possible to define the matrix on the host
only. If the matrix is defined on all nodes, there is no need to
use this function."""
function associate_matrix{Tv <: MUMPSValueDataType, Ti <: MUMPSIntDataType}(mumps :: Mumps{Tv}, A :: SparseMatrixCSC{Tv,Ti})

  n = size(A, 1);
  size(A, 2) == n || throw(MUMPSException("Input matrix must be square"))

  # Symmetric factorization only accesses the lower triangle.
  B = mumps.__sym > 0 ? tril(A) : A;

  # Obtain B in coordinate format.
  nz = nnz(B);
  #   valtype = isreal(B.nzval[1]) ? Float64 : Complex128;
  #   valtype == mumps.valtype || throw(MUMPSException("Inconsistent data type"))
  vals = convert(Array{Tv,1}, B.nzval);       # Necessary?
  irow = convert(Array{Int32,1}, B.rowval);   # Necessary?
  jcol = zeros(Int32, nz, 1);
  for i = 1 : n
    jcol[B.colptr[i] : B.colptr[i+1]-1] = i;
  end

  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    @mumps_call(:mumps_associate_matrix_float, Void,
                (Ptr{Void}, Int32, Int32, Ptr{Float32}, Ptr{Int32}, Ptr{Int32}),
                        id,     n,    nz,         vals,       irow,       jcol);
  elseif Tv == Float64
    @mumps_call(:mumps_associate_matrix_double, Void,
                (Ptr{Void}, Int32, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
                        id,     n,    nz,         vals,       irow,       jcol);
  elseif Tv == Complex64
    @mumps_call(:mumps_associate_matrix_float_complex, Void,
                (Ptr{Void}, Int32, Int32, Ptr{Complex64}, Ptr{Int32}, Ptr{Int32}),
                        id,     n,    nz,           vals,       irow,       jcol);
  else
    @mumps_call(:mumps_associate_matrix_double_complex, Void,
                (Ptr{Void}, Int32, Int32, Ptr{Complex128}, Ptr{Int32}, Ptr{Int32}),
                        id,     n,    nz,            vals,       irow,       jcol);
  end

  mumps.n = n;
  mumps.nnz = mumps.infog[29];
  return mumps;
end

# Associate a generally-typed matrix with a Mumps type. Attempt conversion.
# An InexactError should be raised if, e.g., mumps is Float64 and A is Complex128.
associate_matrix{Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti}) = associate_matrix(mumps, convert(SparseMatrixCSC{Tm,Int64}, A));

# associate_matrix for dense matrices.
associate_matrix{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, A :: Array{Tv,2}) = associate_matrix(mumps, convert(SparseMatrixCSC{Tm,Int64}, sparse(A)));


import Base.LinAlg.factorize

"""Factorize the matrix registered with the `Mumps` instance.
The matrix must have been previously registered with `associate_matrix()`.
After the factorization, the determinant, if requested, is stored in
`mumps.det`. The MUMPS error code is stored in `mumps.err`. """
function factorize{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv})

  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    @mumps_call(:mumps_factorize_float, Void, (Ptr{Void},), id);

    @mumps_call(:mumps_get_info_float, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float32}),
                        id, mumps.infog, mumps.rinfog)
  elseif Tv == Float64
    @mumps_call(:mumps_factorize_double, Void, (Ptr{Void},), id);

    @mumps_call(:mumps_get_info_double, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
                id, mumps.infog, mumps.rinfog)
  elseif Tv == Complex64
    @mumps_call(:mumps_factorize_float_complex, Void, (Ptr{Void},), id);

    @mumps_call(:mumps_get_info_float_complex, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float32}),
                id, mumps.infog, mumps.rinfog)
  else
    @mumps_call(:mumps_factorize_double_complex, Void, (Ptr{Void},), id);

    @mumps_call(:mumps_get_info_double_complex, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
                id, mumps.infog, mumps.rinfog)
  end

  if mumps.icntl[33] == 1
    mumps.det = mumps.rinfog[12] * 2.0^(mumps.infog[34]);
  end
  mumps.err = mumps.infog[1];
  return mumps;
end

"""Register the right-hand side(s) `rhs` with the `Mumps`
object `mumps`. This function makes it possible to define the right-
-hand side(s) on the host only. If the right-hand side(s) are defined
on all nodes, there is no need to use this function."""
function associate_rhs{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}, rhs :: Array{Tv})
  associate_rhs!(mumps, copy(rhs))
end


"""Register the right-hand side(s) `rhs` with the `Mumps`
object `mumps`. This function makes it possible to define the right-
-hand side(s) on the host only. If the right-hand side(s) are defined
on all nodes, there is no need to use this function.

Note that `rhs` will be overwritten with the solution after a call to `solve()`."""
function associate_rhs!{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}, rhs :: Array{Tv})

  n = size(rhs, 1)
  n == mumps.n || throw(MUMPSException("rhs has incompatible dimension"))

  nrhs = size(rhs, 2)

  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    @mumps_call(:mumps_associate_rhs_float, Void,
                (Ptr{Void}, Int32, Ptr{Float32}),
                        id,  nrhs,            x);
  elseif Tv == Float64
    @mumps_call(:mumps_associate_rhs_double, Void,
                (Ptr{Void}, Int32, Ptr{Float64}),
                        id,  nrhs,            x);
  elseif Tv == Complex64
    @mumps_call(:mumps_associate_rhs_float_complex, Void,
                (Ptr{Void}, Int32, Ptr{Complex64}),
                        id,  nrhs,              x);
  else
    @mumps_call(:mumps_associate_rhs_double_complex, Void,
                (Ptr{Void}, Int32, Ptr{Complex128}),
                        id,  nrhs,               x);
  end
  return mumps;
end

# Associate a generally-typed rhs with a Mumps type. Attempt conversion.
# An InexactError should be raised if, e.g., mumps is Float64 and rhs is Complex128.
associate_rhs{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, rhs :: Array{Tv}) = associate_rhs(mumps, convert(Array{Tm}, rhs));


"""Solve the system registered with the `Mumps` object `mumps`.
The matrix and right-hand side(s) must have been previously registered
with `associate_matrix()` and `associate_rhs()`. The optional keyword
argument `transposed` indicates whether the user wants to solve the
forward or transposed system. The solution is stored internally and must
be retrieved with `get_solution()`."""
function solve{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}; transposed :: Bool=false)

  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    @mumps_call(:mumps_solve_float, Void,
                (Ptr{Void}, Int32),
                        id, transposed ? 1 : 0);

    @mumps_call(:mumps_get_info_float, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float32}),
                        id, mumps.infog, mumps.rinfog)
  elseif Tv == Float64
    @mumps_call(:mumps_solve_double, Void,
                (Ptr{Void}, Int32),
                        id, transposed ? 1 : 0);

    @mumps_call(:mumps_get_info_double, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
                        id, mumps.infog, mumps.rinfog)
  elseif Tv == Complex64
    @mumps_call(:mumps_solve_float_complex, Void,
                (Ptr{Void}, Int32),
                        id, transposed ? 1 : 0);

    @mumps_call(:mumps_get_info_float_complex, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float32}),
                        id, mumps.infog, mumps.rinfog)
  else
    @mumps_call(:mumps_solve_double_complex, Void,
                (Ptr{Void}, Int32),
                        id, transposed ? 1 : 0);

    @mumps_call(:mumps_get_info_double_complex, Void,
                (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
                        id, mumps.infog, mumps.rinfog)
  end

  mumps.err = mumps.infog[1];
  return mumps;
end


"""Retrieve the solution of the system solved by `solve()`. This
function makes it possible to ask MUMPS to assemble the final solution
on the host only, and to retrieve it there."""
function get_solution{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv})

  id = reinterpret(Ptr{Void}, mumps.__id)
  if Tv == Float32
    nrhs = int(@mumps_call(:mumps_get_nrhs_float, Int32, (Ptr{Void},), id));

    x = zeros(Float32, mumps.n * nrhs);
    @mumps_call(:mumps_get_solution_float, Void,
                (Ptr{Void}, Ptr{Float32}),
                        id,            x);
  elseif Tv == Float64
    nrhs = int(@mumps_call(:mumps_get_nrhs_double, Int32, (Ptr{Void},), id));

    x = zeros(Float64, mumps.n * nrhs);
    @mumps_call(:mumps_get_solution_double, Void,
                (Ptr{Void}, Ptr{Float64}),
                        id,            x);
  elseif Tv == Complex64
    nrhs = int(@mumps_call(:mumps_get_nrhs_float_complex, Int32, (Ptr{Void},), id));

    x = zeros(Complex64, mumps.n * nrhs);
    @mumps_call(:mumps_get_solution_float_complex, Void,
                (Ptr{Void}, Ptr{Complex64}),
                        id,              x);
  else
    nrhs = int(@mumps_call(:mumps_get_nrhs_double_complex, Int32, (Ptr{Void},), id));

    x = zeros(Complex128, mumps.n * nrhs);
    @mumps_call(:mumps_get_solution_double_complex, Void,
                (Ptr{Void}, Ptr{Complex128}),
                        id,               x);
  end

  return reshape(x, Int(mumps.n), nrhs);
end


# Convenience functions.

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
function factorize{Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti})
  mumps = associate_matrix(mumps, A);  # A will be converted by associate_matrix.
  return factorize(mumps);
end

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
factorize{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, A :: Array{Tv}) = factorize(mumps, convert(SparseMatrixCSC{Tm,Int64}, sparse(A)));


"""Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, rhs :: Array{Tv}; transposed :: Bool=false)
  mumps = associate_rhs(mumps, rhs);  # rhs will be converted by associate_rhs.
  mumps = solve(mumps, transposed=transposed);
  return get_solution(mumps);
end


"""Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve{Tm <: MUMPSValueDataType, Tv <: Number, Tr <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti}, rhs :: Array{Tr}; transposed :: Bool=false)

  mumps = factorize(mumps, A);
  return solve(mumps, rhs, transposed=transposed);
end

solve{Tm <: MUMPSValueDataType, Tv <: Number, Tr <: Number}(mumps :: Mumps{Tm}, A :: Array{Tv,2}, rhs :: Array{Tr}) = solve(mumps, sparse(A), rhs);


"""Combined initialize / analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `sym` indicates the symmetry of `A`.
The solution is retrieved and returned."""
function solve{Tv <: Number, Tr <: Number, Ti <: MUMPSIntDataType}(A :: SparseMatrixCSC{Tv,Ti}, rhs :: Array{Tr}; sym :: Int=mumps_unsymmetric)

  Tm = (Tv <: Complex || Tr <: Complex) ? Complex128 : Float64;  # Could be smarter.
  mumps = Mumps{Tm}(sym, default_icntl, default_cntl64);
  x = solve(mumps, A, rhs);
  finalize(mumps);
  return x;
end

solve{Tv <: Number, Tr <: Number}(A :: Array{Tv,2}, rhs :: Array{Tr}; sym :: Int=mumps_unsymmetric) = solve(sparse(A), rhs, sym=sym);

end  # Module MUMPS
