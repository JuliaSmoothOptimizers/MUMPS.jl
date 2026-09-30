export get_icntl, default_icntl, default_cntl32, default_cntl64

export associate_matrix!, associate_rhs!, set_user_perm!, get_solution, solve!, solve, factorize!

export mumps_unsymmetric, mumps_definite, mumps_symmetric

"""
    associate_matrix!(mumps, n, irow, jcol, vals)
Register the sparse matrix given in coordinate format with the `Mumps` object `mumps`.
This function makes it possible to define the matrix on the host
only. If the matrix is defined on all nodes, there is no need to
use this function.
"""
associate_matrix!(mumps::Mumps, n::Integer, irow::Vector, jcol::Vector, vals::Vector) =
  associate_matrix!(mumps, sparse(n, n, irow, jcol, vals))

"""
    get_solution(mumps) -> x
Retrieve the solution of the system solved by `solve()`. This
function makes it possible to ask MUMPS to assemble the final solution
on the host only, and to retrieve it there.
"""
function get_solution end

get_solution(args...; kwargs...) = get_sol(args...; kwargs...)

"""
    factorize!(mumps)
Factorize the matrix registered with the `Mumps` instance.
The matrix must have been previously registered with `associate_matrix()`.
After the factorization, the determinant, if requested, is stored in
`mumps.det`. The MUMPS error code is stored in `mumps.err`.
"""
function factorize!(mumps::Mumps)
  # suppress_printing!(mumps)
  mumps_factorize!(mumps)
end

"""
    factorize!(mumps,A)
Combined associate_matrix / factorize.
Presume that `A` is available on all nodes.
"""
function factorize!(mumps::Mumps, A::AbstractArray)
  associate_matrix!(mumps, A)
  factorize!(mumps)
  mumps.icntl[33] == 1 && (mumps.det = det(mumps))
  return mumps
end

"""
    solve!(mumps;transposed=false)
Solve the system registered with the `Mumps` object `mumps`.
The matrix and right-hand side(s) must have been previously registered
with `associate_matrix()` and `associate_rhs()`. The optional keyword
argument `transposed` indicates whether the user wants to solve the
forward or transposed system. The solution is stored internally and must
be retrieved with `get_solution()`.
"""
function solve!(mumps::Mumps; transposed::Bool = false)
  # suppress_printing!(mumps)
  transposed && transpose!(mumps)
  mumps_solve!(mumps)
  transposed && transpose!(mumps)
  return mumps
end

"""
    solve(mumps, rhs; transposed=false)
Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned.
"""
function solve(mumps::Mumps, rhs::AbstractArray; transposed::Bool = false)
  associate_rhs!(mumps, rhs)
  solve!(mumps; transposed = transposed)
  return get_sol(mumps)
end

"""
    solve(mumps, A, rhs; transposed=false)
Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned.
"""
function solve(mumps::Mumps, A::AbstractArray, rhs::AbstractArray; transposed::Bool = false)
  # suppress_printing!(mumps)
  factorize!(mumps, A)
  return solve(mumps, rhs; transposed = transposed)
end

"""
    solve(A, rhs; sym=mumps_unsymmetric)
Combined initialize / analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `sym` indicates the symmetry of `A`.
The solution is retrieved and returned.
"""
function solve(
  A::AbstractArray{T},
  rhs::AbstractArray{V};
  sym::Integer = mumps_unsymmetric,
) where {T, V}
  mumps = Mumps{promote_type(T, V)}(sym)
  # suppress_printing!(mumps)
  associate_matrix!(mumps, A)
  associate_rhs!(mumps, rhs)
  solve!(mumps)
  x = get_sol(mumps)
  finalize(mumps)
  return x
end

"""
    finalize(mumps)
Terminate a Mumps instance.
"""
Base.finalize(mumps::Mumps) = finalize!(mumps)

# See Section 6.1 of the MUMPS 5.9.1 User's Guide.
"Default integer parameters."
default_icntl = zeros(Int32, 60);
default_icntl[1] = 6;  # Output stream for error messages
default_icntl[2] = 0;  # Output stream for diagonstics/stats/warnings
default_icntl[3] = 6;  # Output stream for global info on host
default_icntl[4] = 2;  # Output level for errors/warnings/diagnostics
default_icntl[5] = 0;  # 0 = assembled matrix, 1 = elemental format
default_icntl[6] = 7;  # permutation/scaling in analysis (7 = automatic)
default_icntl[7] = 7;  # pivot order for factorization (7 = automatic)
default_icntl[8] = 77;  # scaling in analysis/factorization (77 = automatic)
default_icntl[9] = 1;  # 1: solve Ax=b, otherwise A'x=b
default_icntl[10] = 0;  # max number of iterative refinement steps
default_icntl[11] = 0;  # > 0: return stats collected during solve
default_icntl[12] = 0;  # ordering during analysis
default_icntl[13] = 0;  # 0: use ScaLAPACK (parallel), >0: sequential unless workers > value
default_icntl[14] = 20;  # % workspace increase during analysis/fact
default_icntl[15] = 0;  # analysis by blocks: 0 = off, 1 = user (BLKPTR), <0 = block size -value
default_icntl[16] = 0;  # number of OpenMP threads set by MUMPS (0 = do not change)
default_icntl[17] = 0;  # use some MPI processes as OpenMP resources (0 = off, MUMPS ≥ 5.9)
default_icntl[18] = 0;  # 0 = matrix assembled on host
default_icntl[19] = 0;  # 1 = return Schur complement on host
default_icntl[20] = 0;  # 0 = dense rhs, 1-3 = sparse rhs, 10-11 = distributed rhs
default_icntl[21] = 0;  # 0 = solution overwrites rhs, 1 = keep distributed
default_icntl[22] = 0;  # 0 = in core, 1 = out of core
default_icntl[23] = 0;  # max working memory
default_icntl[24] = 0;  # 0: null pivot=error
default_icntl[25] = 0;  # -1: compute nullspace basis
default_icntl[26] = 0;  # condense rhs on Schur variables (see 19)
default_icntl[27] = -32;  # blocking size for multiple rhs (<0: automatic)
default_icntl[28] = 1;  # 1: sequential analysis, 2: parallel, 0: automatic
default_icntl[29] = 0;  # ordering for parallel analysis (see 28)
default_icntl[30] = 0;  # compute entries of the inverse
default_icntl[31] = 0;  # discard factors after factorization (1: all, 2: L only for unsymmetric)
default_icntl[32] = 0;  # 1: forward elimination during factorization
default_icntl[33] = 0;  # compute determinant
default_icntl[34] = 0;  # save/restore: 1 = keep OOC files when deleting saved data
default_icntl[35] = 0;  # BLR: 0 = off, 1 = auto, 2 = facto + solve, 3 = facto only
default_icntl[36] = 1;  # BLR variant: 0 = UFSC, 1 = UCFS (default since MUMPS 5.9)
default_icntl[37] = 0;  # BLR compression of contribution blocks
default_icntl[38] = 600;  # estimated compression rate of LU factors (per mille)
default_icntl[39] = 500;  # estimated compression rate of contribution blocks (per mille)
default_icntl[40] = 0;  # mixed/adaptive precision BLR (0 = off, MUMPS ≥ 5.9)
# default_icntl[41:46] are not used.
default_icntl[47] = 0;  # 1: single precision facto in double precision instance (≥ 5.9)
default_icntl[48] = 1;  # 1: multithreading with tree parallelism (L0-threads)
default_icntl[49] = 0;  # compact workarray S at the end of factorization (0 = off, 1 or 2 = on)
# default_icntl[50] is not used.
default_icntl[51] = 0;  # offload activities to GPUs (0 = off, MUMPS ≥ 5.9)
# default_icntl[52:55] are not used.
default_icntl[56] = 0;  # 1: rank-revealing factorization (null space detection)
# default_icntl[57] is not used.
default_icntl[58] = 2;  # symbolic factorization: 1 = quotient graph, 2 = column counts
# default_icntl[59:60] are not used.

# See Section 6.2 of the MUMPS 5.9.1 User's Guide.
# cntl[1] = -1 lets MUMPS choose the pivoting threshold automatically (depends on SYM and on
# rank-revealing, ICNTL(56)). cntl[7] is the BLR dropping parameter (0 = no compression).
"Default single precision real parameters"
default_cntl32 = zeros(Float32, 15);
default_cntl32[1] = -1;    # relative threshold for numerical pivoting
default_cntl32[2] = sqrt(eps(Float32));  # tolerance for iterative refinement
default_cntl32[3] = 0.0;  # threshold to detect null pivots
default_cntl32[4] = -1.0;  # threshold for static pivoting (<0: disable)
default_cntl32[5] = 0.0;  # what null pivots are reset to
default_cntl32[7] = 0.0;  # BLR dropping parameter ε (0 = full precision)
# default_cntl32[6] and default_cntl32[8-15] are not used.

"Default double precision real parameters"
default_cntl64 = zeros(Float64, 15);
default_cntl64[1] = -1;    # relative threshold for numerical pivoting
default_cntl64[2] = sqrt(eps(Float64));  # tolerance for iterative refinement
default_cntl64[3] = 0.0;  # threshold to detect null pivots
default_cntl64[4] = -1.0;  # threshold for static pivoting (<0: disable)
default_cntl64[5] = 0.0;  # what null pivots are reset to
default_cntl64[7] = 0.0;  # BLR dropping parameter ε (0 = full precision)
# default_cntl64[6] and default_cntl64[8-15] are not used.

"""
    get_icntl(; det=false, verbose=false, ooc=false, itref=0, user_perm=false,
              blr=0, nthreads=0, tree_parallelism=true, rank_revealing=false,
              single_precision_factorization=false, null_pivots=false)

Obtain an array of integer control parameters.

* `det`: compute the determinant (ICNTL(33));
* `verbose`: print MUMPS output (ICNTL(1:4));
* `ooc`: store factors out of core (ICNTL(22));
* `itref`: maximum number of iterative refinement steps (ICNTL(10));
* `user_perm`: use a user-supplied permutation, see [`set_user_perm!`](@ref) (ICNTL(7));
* `blr`: Block Low-Rank factorization, 0 = off, 1 = automatic, 2 = factorization and solve,
  3 = factorization only (ICNTL(35)); the low-rank tolerance is CNTL(7);
* `nthreads`: number of OpenMP threads MUMPS should use, 0 = leave unchanged (ICNTL(16));
* `tree_parallelism`: multithreading exploiting tree parallelism, a.k.a. L0-threads (ICNTL(48));
* `rank_revealing`: rank-revealing factorization for null space detection (ICNTL(56));
* `single_precision_factorization`: factorize in single precision within a `Float64` or
  `ComplexF64` instance (ICNTL(47), MUMPS ≥ 5.9); combine with `itref` to recover double
  precision accuracy;
* `null_pivots`: detect null pivot rows (ICNTL(24)).
"""
function get_icntl(;
  det::Bool = false,       # Compute determinant.
  verbose::Bool = false,   # Output intermediate info.
  ooc::Bool = false,       # Store factors out of core.
  itref::Int = 0,          # Max steps of iterative refinement.
  user_perm::Bool = false, # Use user-supplied permutation.
  blr::Int = 0,            # Block Low-Rank factorization.
  nthreads::Int = 0,       # Number of OpenMP threads (0 = unchanged).
  tree_parallelism::Bool = true, # L0-threads multithreading.
  rank_revealing::Bool = false,  # Rank-revealing factorization.
  single_precision_factorization::Bool = false, # Mixed precision factorization.
  null_pivots::Bool = false, # Null pivot row detection.
)
  0 ≤ blr ≤ 3 || throw(ArgumentError("blr must be 0, 1, 2 or 3, got $blr"))
  nthreads ≥ 0 || throw(ArgumentError("nthreads must be nonnegative, got $nthreads"))
  icntl = default_icntl[:]
  icntl[33] = det ? 1 : 0
  if !verbose
    icntl[1:4] .= 0
  end
  icntl[22] = ooc ? 1 : 0
  icntl[10] = itref
  icntl[7] = user_perm ? 1 : 7  # 1 = user-supplied, 7 = automatic
  icntl[35] = blr
  icntl[16] = nthreads
  icntl[48] = tree_parallelism ? 1 : 0
  icntl[56] = rank_revealing ? 1 : 0
  icntl[47] = single_precision_factorization ? 1 : 0
  icntl[24] = null_pivots ? 1 : 0
  return icntl
end

# Symbols for symmetry
const mumps_unsymmetric = 0
const mumps_definite = 1
const mumps_symmetric = 2

"""
    const mumps_unsymmetric
Constant indicating that a general unsymmetric matrix will be
analyzed and factorized
"""
mumps_unsymmetric

"""
    const mumps_definite
Constant indicating that a symmetric definite matrix will be
analyzed and factorized
"""
mumps_definite

"""
    const mumps_symmetric
Constant indicating that a general symmetric matrix will be
analyzed and factorized
"""
mumps_symmetric
