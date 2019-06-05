export get_icntl,
default_icntl,
default_cntl32,
default_cntl64

export associate_matrix!,
associate_rhs!,
get_solution,
solve!,
solve,
factorize!

export mumps_unsymmetric,
mumps_definite,
mumps_symmetric


"""
    associate_matrix!(mumps,A)
Register the matrix `A` with the `Mumps` object `mumps`.
This function makes it possible to define the matrix on the host
only. If the matrix is defined on all nodes, there is no need to
use this function.
"""
associate_matrix!(args...;kwargs...) = provide_matrix!(args...;kwargs...)
associate_matrix!(mumps::Mumps, n::Integer, irow::Vector, jcol::Vector, vals::Vector) =
    associate_matrix!(mumps,sparse(n,n,irow,jcol,vals))


"""
    associate_rhs!(mumps,rhs)
Register the right-hand side(s) `rhs` with the `Mumps`
object `mumps`. This function makes it possible to define the right-
-hand side(s) on the host only. If the right-hand side(s) are defined
on all nodes, there is no need to use this function.
"""
associate_rhs! = provide_rhs!

"""
    get_solution(mumps) -> x
Retrieve the solution of the system solved by `solve()`. This
function makes it possible to ask MUMPS to assemble the final solution
on the host only, and to retrieve it there.
"""
get_solution = get_sol

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
function factorize!(mumps::Mumps,A::AbstractArray)
    provide_matrix!(mumps,A)
    factorize!(mumps)
    mumps.icntl[33]==1 ? mumps.det = det(mumps) : nothing
    return nothing
end


"""
    solve!(mumps;transposed=false)
Solve the system registered with the `Mumps` object `mumps`.
The matrix and right-hand side(s) must have been previously registered
with `associate_matrix()` and `associate_rhs()`. The optional keyword
argument `transposed` indicates whether the user wants to solve the
forward or transposed system. The solution is stored internally and must
be retrieved with `get_solution()`."""
function solve!(mumps::Mumps; transposed::Bool=false)
    # suppress_printing!(mumps)
    transposed ? transpose!(mumps) : nothing
    mumps_solve!(mumps)
    transposed ? transpose!(mumps) : nothing
    return nothing
end


"""
    solve(mumps, rhs; transposed=false)
Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve(mumps::Mumps,rhs::AbstractArray; transposed::Bool=false)
    provide_rhs!(mumps, rhs)
    solve!(mumps; transposed=transposed)
    return get_sol(mumps)
end
"""
    solve(mumps,A,rhs; transposed=false)
Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve(mumps::Mumps, A::AbstractArray, rhs::AbstractArray; transposed::Bool=false)
    # suppress_printing!(mumps)
    factorize!(mumps,A)
    return solve(mumps, rhs; transposed=transposed)
end
"""
    solve(A,rhs;sym=mumps_unsymmetric)
Combined initialize / analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `sym` indicates the symmetry of `A`.
The solution is retrieved and returned.
"""
function solve(A::AbstractArray{T}, rhs::AbstractArray{V}; sym::Integer=mumps_unsymmetric) where {T,V}
    mumps = Mumps{promote_type(T,V)}(sym)
    # suppress_printing!(mumps)
    provide_matrix!(mumps,A)
    provide_rhs!(mumps,rhs)
    solve!(mumps)
    return get_sol(mumps)
end


"""
    finalize(mumps)
Terminate a Mumps instance.
"""
Base.finalize(mumps::Mumps) = finalize!(mumps)


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


"""
    get_icntl(;det=false,verbose=false,ooc=false,itref=0)
Obtain an array of integer control parameters.
"""
function get_icntl(;
                   det :: Bool=false,       # Compute determinant.
                   verbose :: Bool=false,   # Output intermediate info.
                   ooc :: Bool=false,       # Store factors out of core.
                   itref :: Int=0,          # Max steps of iterative refinement.
                   )
  icntl = default_icntl[:]
  icntl[33] = det ? 1 : 0
  if !verbose
    icntl[1:4] .= 0
  end
  icntl[22] = ooc ? 1 : 0
  icntl[10] = itref
  return icntl
end


# Symbols for symmetry
const mumps_unsymmetric = 0
const mumps_definite    = 1
const mumps_symmetric   = 2


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
