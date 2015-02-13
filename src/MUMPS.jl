module MUMPS

export default_icntl, default_cntl, Mumps, finalize, factorize, solve,
       associate_matrix, associate_rhs, get_solution,
       mumps_unsymmetric, mumps_definite, mumps_symmetric,
       MUMPSException

using MPI
using Docile
@docstrings(manual = ["../doc/manual.md"])

# libjmumps.dylib should be on your LD_LIBRARY_PATH.
mumps_lib = "libjmumps";
macro mumps_call(func, args...)
  quote
    ccall(($func, $mumps_lib), $(args...))
  end
end


@doc "Exception type raised in case of error." ->
type MUMPSException <: Exception
  msg :: ASCIIString
end

MUMPSDataType = Union(Type{Float64}, Type{Complex128})


# See MUMPS User's Manual Section 5.1.
@doc "Default integer parameters." ->
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
@doc "Default real parameters" ->
default_cntl = zeros(Float64, 15);
default_cntl[1] = -1;    # relative threshold for numerical pivoting
default_cntl[2] = sqrt(eps(Float64));  # tolerance for iterative refinement
default_cntl[3] =  0.0;  # threshold to detect null pivots
default_cntl[4] = -1.0;  # threshold for static pivoting (<0: disable)
default_cntl[5] =  0.0;  # what null pivots are reset to
# default_cntl[6-15] are not used.

# Symbols for symmetry
@doc """Constant indicating that a general unsymmetric matrix will be
analyzed and factorized""" ->
mumps_unsymmetric = 0;

@doc """Constant indicating that a symmetric definite matrix will be
analyzed and factorized""" ->
mumps_definite    = 1;

@doc """Constant indicating that a general symmetric matrix will be
analyzed and factorized""" ->
mumps_symmetric   = 2;


@doc """Abstract type representing a factorization with MUMPS.
All constructor arguments are optional. By default a general
unsymmetric matrix will be analyzed/factorized with default
integer and real parameters""" ->
type Mumps
  __id    :: Ptr{Void}         # Pointer to MUMPS struct. Do not touch.
  __sym   :: Int32             # Value of sym used by Mumps.
  icntl   :: Array{Int32,1}    # Integer control parameters.
  cntl    :: Array{Float64,1}  # Real control parameters.
  n       :: Int32             # Order of the matrix factorized.
  valtype :: MUMPSDataType     # Real or complex.
  infog   :: Array{Int32,1}
  rinfog  :: Array{Float64,1}
  nnz     :: Int               # Number of nonzeros in factors.
  det     :: Int
  err     :: Int

  # Constructor.
  function Mumps(sym   :: Int=mumps_unsymmetric;
                 icntl :: Array{Int32,1}=default_icntl,
                 cntl  :: Array{Float64,1}=default_cntl)
    # sym = 0 (unsymmetric), 1 (symmetric definite), 2 (symmetric indefinite).
    # There is currently no facility to exploit hermicity of complex matrices.

    MPI.Initialized() || throw(MUMPSException("Initialize MPI first"));

    # Set default pivot threshold if required.
    if cntl[1] == -1
      cntl[1] = (sym == mumps_definite) ? 0.0 : 0.01
    end

    id = @mumps_call(:mumps_initialize_double, Ptr{Void},
                     (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);

    id == C_NULL && throw(MUMPSException("Error allocating MUMPS structure"))

    infog = zeros(Int32, 40);
    rinfog = zeros(Float64, 20);

    self = new(id, int32(sym), icntl, cntl, 0, Float64, infog, rinfog, 0, 0, 0);
    finalizer(self, finalize);  # Destructor.
    return self;
  end
end


# Convenience constructor.
# Note that every argument is a keyword.
@doc """Convenience constructor. All arguments are optional.
* sym: one of mumps_unsymmetric, mumps_definite, mumps_symmetric
* verbose: true or false
* det: true or false, compute determinant
* ooc: true or false, store factors out of core
* itref: number of iterative refinement steps
* cntl: real parameters array
""" ->
function Mumps(;sym :: Int=mumps_unsymmetric,
               det :: Bool=false,       # Compute determinant.
               verbose :: Bool=false,   # Output intermediate info.
               ooc :: Bool=false,       # Store factors out of core.
               itref :: Int=0,          # Max steps of iterative refinement.
               cntl :: Array{Float64,1}=default_cntl
               )

  icntl = default_icntl[:];
  icntl[33] = det ? 1 : 0;
  if !verbose
    icntl[1:4] = 0;
  end
  icntl[22] = ooc ? 1 : 0;
  icntl[10] = itref;

  return Mumps(sym, icntl=icntl, cntl=cntl);
end


@doc "Terminate a Mumps instance." ->
function finalize(mumps :: Mumps)
  @mumps_call(:mumps_finalize_double, Void, (Ptr{Void},), mumps.__id);
  mumps.__id = C_NULL;
end


@doc """Register the matrix `A` with the `Mumps` object `mumps`.
This function makes it possible to define the matrix on the host
only. If the matrix is defined on all nodes, there is no need to
use this function.""" ->
function associate_matrix(mumps :: Mumps, A :: SparseMatrixCSC)

  n = size(A, 1);
  size(A, 2) == n || throw(MUMPSException("Input matrix must be square"))

  # Symmetric factorization only accesses the lower triangle.
  B = mumps.__sym > 0 ? tril(A) : A;

  # Obtain B in coordinate format.
  nz = nnz(B);
  valtype = isreal(B.nzval[1]) ? Float64 : Complex128;
  vals = convert(Array{valtype,1}, B.nzval);
  irow = convert(Array{Int32,1}, B.rowval);
  jcol = zeros(Int32, nz, 1);
  for i = 1 : n
    jcol[B.colptr[i] : B.colptr[i+1]-1] = i;
  end

  @mumps_call(:mumps_associate_matrix_double, Void,
              (Ptr{Void}, Int32, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
              mumps.__id,     n,    nz,         vals,       irow,       jcol);

  mumps.n = n;
  mumps.valtype = valtype;
  mumps.nnz = mumps.infog[29];
  return;
end


associate_matrix(mumps :: Mumps, A :: Array{Float64}) = associate_matrix(mumps, sparse(A));


import Base.LinAlg.factorize

@doc """Factorize the matrix registered with the `Mumps` instance.
The matrix must have been previously registered with `associate_matrix()`.
After the factorization, the determinant, if requested, is stored in
`mumps.det`. The MUMPS error code is stored in `mumps.err`. """ ->
function factorize(mumps :: Mumps)

  @mumps_call(:mumps_factorize_double, Void, (Ptr{Void},), mumps.__id);

  @mumps_call(:mumps_get_info_double, Void,
              (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
              mumps.__id, mumps.infog, mumps.rinfog)

  if mumps.icntl[33] == 1
    mumps.det = mumps.rinfog[12] * 2^(mumps.infog[34]);
  end
  mumps.err = mumps.infog[1];
  return;
end


@doc """Register the right-hand side(s) `rhs` with the `Mumps`
object `mumps`. This function makes it possible to define the right-
-hand side(s) on the host only. If the right-hand side(s) are defined
on all nodes, there is no need to use this function.""" ->
function associate_rhs(mumps :: Mumps, rhs :: Array{Float64})

  n = size(rhs, 1);
  n == mumps.n || throw(MUMPSException("rhs has incompatible dimension"))
  typeof(rhs[1]) == mumps.valtype || throw(MUMPSException("rhs has incompatible data type"))

  nrhs = size(rhs, 2);
  x = rhs[:];  # Make a copy; will be overwritten with solution.

  @mumps_call(:mumps_associate_rhs_double, Void,
              (Ptr{Void}, Int32, Ptr{Float64}),
              mumps.__id,  nrhs,            x);
  return;
end


@doc """Solve the system registered with the `Mumps` object `mumps`.
The matrix and right-hand side(s) must have been previously registered
with `associate_matrix()` and `associate_rhs()`. The optional keyword
argument `transposed` indicates whether the user wants to solve the
forward or transposed system. The solution is stored internally and must
be retrieved with `get_solution()`.""" ->
function solve(mumps :: Mumps; transposed :: Bool=false)

  @mumps_call(:mumps_solve_double, Void,
              (Ptr{Void}, Int32),
              mumps.__id, transposed ? 1 : 0);

  @mumps_call(:mumps_get_info_double, Void,
              (Ptr{Void}, Ptr{Int32},  Ptr{Float64}),
              mumps.__id, mumps.infog, mumps.rinfog)

  mumps.err = mumps.infog[1];
  return;
end


@doc """Retrieve the solution of the system solved by `solve()`. This
function makes it possible to ask MUMPS to assemble the final solution
on the host only, and to retrieve it there.""" ->
function get_solution(mumps :: Mumps)
  nrhs = int(@mumps_call(:mumps_get_nrhs_double, Int32, (Ptr{Void},), mumps.__id));

  x = zeros(Float64, mumps.n * nrhs);
  @mumps_call(:mumps_get_solution_double, Void,
              (Ptr{Void}, Ptr{Float64}),
              mumps.__id,            x);

  x = zeros(Float64, mumps.n * nrhs);
  @mumps_call(:mumps_get_solution, Void,
              (Ptr{Void}, Ptr{Float64}),
              mumps.__id,            x);
  return reshape(x, int(mumps.n), nrhs);
end


# Convenience functions.

@doc """Combined associate_matrix / factorize.
Presume that `A` is available on all nodes.""" ->
function factorize(mumps :: Mumps, A :: SparseMatrixCSC)
  associate_matrix(mumps, A);
  factorize(mumps);
  return;
end

@doc """Combined associate_matrix / factorize.
Presume that `A` is available on all nodes.""" ->
factorize(mumps :: Mumps, A :: Array{Float64}) = factorize(mumps, sparse(A));


@doc meta("""Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned.""", returns=(Array{Float64},)) ->
function solve(mumps :: Mumps, rhs :: Array{Float64}; transposed :: Bool=false)
  associate_rhs(mumps, rhs);
  solve(mumps, transposed=transposed);
  return get_solution(mumps);
end


@doc meta("""Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned.""", returns=(Array{Float64},)) ->
function solve(mumps :: Mumps, A :: SparseMatrixCSC, rhs :: Array{Float64};
  transposed :: Bool=false)

  factorize(mumps, A);
  return solve(mumps, rhs, transposed=transposed);
end

solve(mumps :: Mumps, A :: Array{Float64}, rhs :: Array{Float64}) = solve(mumps, sparse(A), rhs);


@doc meta("""Combined initialize / analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `sym` indicates the symmetry of `A`.
The solution is retrieved and returned.""", returns=(Array{Float64},)) ->
function solve(A :: SparseMatrixCSC, rhs :: Array{Float64};
               sym :: Int=mumps_unsymmetric)

  mumps = Mumps(sym);
  x = solve(mumps, A, rhs);
  finalize(mumps);
  return x;
end

solve(A :: Array{Float64}, rhs :: Array{Float64}; sym :: Int=mumps_unsymmetric) = solve(sparse(A), rhs, sym=sym);

end  # Module MUMPS
