# this file provides the low-level interface with the MUMPS library
# by controlling access to the pointers contained in the Mumps object.
# Many functions are unsafe and are marked as such.
# None of these functions change the JOB parameter.

"""
    invoke_mumps_unsafe!(mumps)

Call the appropriate mumps C-library, passing to it the Mumps object `mumps`.

This is a low-level function, meaning that you have complete control over what
operations are done, based on the MUMPS manual.

Be warned, a direct call can crash Julia if `mumps` is not appropriately
initialized.

See also: [`invoke_mumps!`](@ref)
"""
invoke_mumps_unsafe!

for (fname, lname, elty, subty) in (("smumps_c", libsmumps, Float32   , Float32),
                                    ("dmumps_c", libdmumps, Float64   , Float64),
                                    ("cmumps_c", libcmumps, ComplexF32, Float32),
                                    ("zmumps_c", libzmumps, ComplexF64, Float64))

  @eval begin
    function invoke_mumps_unsafe!(mumps :: Mumps{$elty,$subty})
      MPI.Initialized() || throw(MUMPSException("must call MPI.Init() exactly once before calling mumps"))
      ccall(($fname, $lname), Cvoid, (Ref{Mumps{$elty,$subty}},), mumps)
      mumps.err = mumps.infog[1]
      return mumps
    end
  end
end

"""
    invoke_mumps!(mumps)

Call the appropriate mumps C-library, passing to it the Mumps object `mumps`,
but checking to make sure `mumps` has been initialized first, so that it's safe.

This is a low-level function, meaning that you have complete control over what
operations are done, based on the MUMPS manual.

Be warned, a direct call can crash Julia if `mumps` is not appropriately
initialized.

See also: [`invoke_mumps_unsafe!`](@ref)
"""
invoke_mumps!(mumps::Mumps) = is_finalized(mumps) ? mumps : invoke_mumps_unsafe!(mumps)
# function check_finalized(mumps::Mumps)
# if mumps._finalized
# throw(MUMPSException("Mumps object already finalized"))
# end
# end

"""
    set_icntl!(mumps,i,val; [displaylevel=1])

Set the integer control parameters according to ICNTL[i]=val

See also: [`display_icntl`](@ref)
"""
function set_icntl!(mumps::Mumps, i::Integer, val::Integer; displaylevel = mumps.icntl[4] - 1)
  icntl = mumps.icntl
  mumps.icntl = (icntl[1:(i - 1)]..., convert(MUMPS_INT, val), icntl[(i + 1):end]...)
  displaylevel > 0 && display_icntl(stdout, mumps.icntl, i, val)
  return mumps
end

"""
    set_cntl!(mumps,i,val; [displaylevel=1])

Set the real/complex control parameters according to CNTL[i]=val

See also: [`display_cntl`](@ref)
"""
function set_cntl!(
  mumps::Mumps{TC, TR},
  i::Integer,
  val::AbstractFloat;
  displaylevel = mumps.icntl[4] - 1,
) where {TC, TR}
  cntl = mumps.cntl
  mumps.cntl = (cntl[1:(i - 1)]..., convert(TR, val), cntl[(i + 1):end]...)
  displaylevel > 0 && display_cntl(stdout, mumps.cntl, i, val)
  return mumps
end

"""
    set_job!(mumps,job)

Set the phase to `job`. See MUMPS manual for options.
"""
function set_job!(mumps::Mumps, i)
  mumps.job = i
  return mumps
end

"""
    set_save_dir!(mumps,dir)

set name of directory in which to store out-of-core files.
"""
function set_save_dir!(mumps, dir::String)
  length(dir) ≤ 255 ||
    throw(MUMPSException("directory name has $(length(dir)) characters, must be ≤ 255"))
  i = length(dir + 1)
  save_dir = mumps.save_dir
  mumps.save_dir = (dir..., '\0', save_dir[(i + 2):end]...)
  return mumps
end

"""
    set_save_prefix!(mumps,dir)

prefix for out-of-core files.
"""
function set_save_prefix!(mumps, prefix::String)
  length(prefix) ≤ 255 ||
    throw(MUMPSException("prefix name has $(length(prefix)) characters, must be ≤ 255"))
  i = length(prefix + 1)
  save_prefix = mumps.save_prefix
  mumps.save_prefix = (prefix..., '\0', save_prefix[(i + 2):end]...)
  return mumps
end

"""
    associate_matrix!(mumps,A)

Register the square matrix `A` to a `mumps` object. It internally converts
`A` to be consistent with the ICNTL[5] setting.

If needed, it tries to convert element type of `A` to be consistent with
type of `mumps`, throwing a warning in this case.

See also: [`associate_rhs!`](@ref)
"""
function associate_matrix!(mumps::Mumps{T}, A::AbstractArray{TA}) where {T, TA}
  size(A, 1) == size(A, 2) ||
    throw(MUMPSException("input matrix must be square, but it is $(size(A,1))×$(size(A,2))"))
  T == TA || (@warn "matrix with element type $TA: will attempt conversion to Mumps type $T")
  if is_matrix_assembled(mumps)
    typeof(A) <: SparseMatrixCSC ||
      (@warn "matrix is dense, but ICNTL[5]=$(mumps.icntl[5]) indicates assembled. attempting to convert matrix to sparse")
    return _associate_matrix_assembled!(mumps, convert(SparseMatrixCSC, A))
  else
    typeof(A) <: SparseMatrixCSC &&
      (@warn "matrix is sparse, but ICNTL[5]=$(mumps.icntl[5]) indicates elemental. attempting to convert matrix to dense")
    return _associate_matrix_elemental!(mumps, convert(Matrix, A))
  end
end
function _associate_matrix_assembled!(mumps::Mumps, A::SparseMatrixCSC{T}) where {T}
  if is_matrix_distributed(mumps)
    return _associate_matrix_assembled_distributed!(mumps, A)
  end
  return _associate_matrix_assembled_centralized!(mumps, A)
end
function _associate_matrix_assembled_centralized!(mumps::Mumps{T}, A::SparseMatrixCSC) where {T}
  if is_symmetric(mumps)
    I, J, V = findnz(triu(A))
  else
    I, J, V = findnz(A)
  end
  irn, jcn, a = convert.((Array{MUMPS_INT}, Array{MUMPS_INT}, Array{T}), (I, J, V))
  mumps.irn, mumps.jcn, mumps.a = pointer.((irn, jcn, a))
  mumps.n = A.n
  mumps.nnz = length(V)
  append!(mumps._gc_haven, [Ref(irn), Ref(jcn), Ref(a)])
  return mumps
end
function _associate_matrix_assembled_distributed!(mumps::Mumps{T}, A::SparseMatrixCSC) where {T}
  throw(MUMPSException("not written yet."))
end
function _associate_matrix_elemental!(mumps::Mumps{T}, A::Array) where {T}
  mumps.n = size(A, 1)
  mumps.nelt = 1
  eltptr = convert.(MUMPS_INT, [1, mumps.n + 1])
  eltvar = convert.(MUMPS_INT, collect(1:(mumps.n)))
  mumps.eltptr, mumps.eltvar = pointer.(eltptr, eltvar)
  if is_symmetric(mumps)
    a_elt = convert.(T, [A[i, j] for i ∈ 1:(mumps.n) for j ∈ 1:i])
  else
    a_elt = convert.(T, A[:])
  end
  mumps.a_elt = pointer(a_elt)
  append!(mumps._gc_haven, [Ref(eltptr), Ref(eltvar), Ref(a_elt)])
  return mumps
end

"""
    associate_rhs!(mumps, rhs)

Register a dense or sparse RHS matrix or vector `rhs` to a `mumps` object. It internally converts
`rhs` to be consistent with the ICNTL[20] setting, and additionally allocates
`mumps.rhs` according to the ICNTL[21] setting.

If needed, it tries to convert element type of `rhs` to be consistent with
type of `mumps`, throwing a warning in this case.

See also: [`associate_matrix!`](@ref)
"""
function associate_rhs!(mumps::Mumps, rhs::AbstractMatrix)
  if is_rhs_dense(mumps)
    return associate_rhs_dense!(mumps, rhs)
  end
  return associate_rhs_sparse!(mumps, rhs)
end
function associate_rhs_sparse!(mumps::Mumps{T}, rhs::AbstractMatrix) where {T}
  rhs = convert(SparseMatrixCSC, rhs)

  mumps.nz_rhs = length(rhs.nzval)
  mumps.nrhs = size(rhs, 2)

  rhs_sparse = convert.(T, rhs.nzval)
  irhs_sparse = convert.(MUMPS_INT, rhs.rowval)
  irhs_ptr = convert.(MUMPS_INT, rhs.colptr)

  mumps.rhs_sparse = pointer(rhs_sparse)
  mumps.irhs_sparse = pointer(irhs_sparse)
  mumps.irhs_ptr = pointer(irhs_ptr)
  if is_sol_central(mumps)
    y = fill(convert(T, NaN), prod(size(rhs)))
    mumps.rhs = pointer(y)
    mumps.lrhs = size(rhs, 1)
    push!(mumps._gc_haven, Ref(y))
  end
  append!(mumps._gc_haven, [Ref(rhs_sparse), Ref(irhs_sparse), Ref(irhs_ptr)])
  return mumps
end
function associate_rhs_dense!(mumps::Mumps{T}, rhs::AbstractMatrix) where {T}
  y = convert(Matrix{T}, rhs)[:]
  mumps.rhs = pointer(y)
  mumps.lrhs = size(rhs, 1)
  mumps.nrhs = size(rhs, 2)
  push!(mumps._gc_haven, Ref(y))
  return mumps
end
associate_rhs!(mumps::Mumps, rhs::AbstractVector) = associate_rhs!(mumps, repeat(rhs, 1, 1))

"""
    get_rhs!(x,mumps)

Retrieve right hand side from `mumps`, storing it in pre-allocated `x`

See also: [`get_rhs`](@ref), [`get_sol!`](@ref), [`get_sol`](@ref)
"""
function get_rhs!(x, mumps::Mumps)
  if !is_finalized(mumps)
    if has_rhs(mumps)
      get_rhs_unsafe!(x, mumps)
    else
      @warn "no rhs registered"
    end
  end
  return x
end
function get_rhs_unsafe!(x::SparseMatrixCSC, mumps::Mumps)
  is_rhs_dense(mumps) &&
    throw(MUMPSException("rhs is dense, target is sparse. try with dense target"))
  for i ∈ LinearIndices(x.colptr)
    x.colptr[i] = unsafe_load(mumps.irhs_ptr, i)
  end
  for i ∈ LinearIndices(x.rowval)
    x.rowval[i] = unsafe_load(mumps.irhs_sparse, i)
    x.nzval[i] = unsafe_load(mumps.rhs_sparse, i)
  end
  return x
end
function get_rhs_unsafe!(x::Union{SubArray, Array}, mumps::Mumps)
  is_rhs_dense(mumps) ||
    throw(MUMPSException("rhs is sparse, target is dense. try with sparse target"))
  for i ∈ LinearIndices(x)
    x[i] = unsafe_load(mumps.rhs, i)
  end
  return x
end
"""
    get_rhs(mumps) -> y

Retrieve right hand side from `mumps`

See also: [`get_rhs!`](@ref), [`get_sol!`](@ref), [`get_sol`](@ref)
"""
function get_rhs(mumps::Mumps{T}) where {T}
  n = mumps.nrhs
  if !is_rhs_dense(mumps)
    m = mumps.n
    colptr = ones(MUMPS_INT, mumps.nrhs + 1)
    rowval = ones(MUMPS_INT, mumps.nz_rhs)
    nzval = Array{T}(undef, mumps.nz_rhs)
    x = SparseMatrixCSC(m, n, colptr, rowval, nzval)
  else
    m = mumps.lrhs
    x = Array{T}(undef, m, n)
  end
  get_rhs!(x, mumps)
end

"""
    get_sol!(x,mumps)

Retrieve solution `x` from `mumps` into pre-allocated array.

See also: [`get_rhs!`](@ref), [`get_rhs`](@ref), [`get_sol`](@ref)
"""
function get_sol!(x::Union{SubArray, Array}, mumps::Mumps)
  if !is_finalized(mumps)
    if mumps.job ∉ [3, 5, 6]
      @error "mumps has not passed through a solution phase"
    end
    if has_rhs(mumps)
      get_sol_unsafe!(x, mumps)
    else
      @warn "no rhs registered"
    end
  end
  return x
end
function get_sol_unsafe!(x::Union{SubArray, Array}, mumps::Mumps)
  for i ∈ LinearIndices(x)
    x[i] = unsafe_load(mumps.rhs, i)
  end
  return x
end
"""
    get_sol(mumps) -> x

Retrieve solution from `mumps`

See also: [`get_rhs!`](@ref), [`get_rhs`](@ref), [`get_sol!`](@ref)
"""
function get_sol(mumps::Mumps{T}) where {T}
  x = Array{T}(undef, mumps.lrhs, mumps.nrhs)
  get_sol!(x, mumps)
end

"""
    get_schur_complement!(S,mumps)

Retrieve Schur complement matrix from `mumps` into pre-allocated `S`

See also: [`get_schur_complement`](@ref), [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function get_schur_complement!(S, mumps::Mumps)
  has_schur(mumps) || throw(MUMPSException("schur complement not yet allocated."))
  get_schur_complement_unsafe!(S, mumps)
end
function get_schur_complement_unsafe!(S, mumps::Mumps)
  for i ∈ LinearIndices(S)
    S[i] = unsafe_load(mumps.schur, i)
  end
  return S
end
"""
    get_schur_complement(mumps) -> S

Retrieve Schur complement matrix `S` from `mumps`

See also: [`get_schur_complement!`](@ref), [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function get_schur_complement(mumps::Mumps{T}) where {T}
  S = Array{T}(undef, mumps.size_schur, mumps.size_schur)
  get_schur_complement!(S, mumps)
end

"""
    set_schur_centralized_by_column!(mumps,schur_inds)

Set up Schur complement matrix calculation for the "centralized by column"
method suggested in the MUMPS manual

See also: [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function set_schur_centralized_by_column!(mumps::Mumps{T}, schur_inds::AbstractArray{Int}) where {T}
  mumps.size_schur = length(schur_inds)
  listvar_schur = convert.(MUMPS_INT, schur_inds)
  mumps.listvar_schur = pointer(listvar_schur)
  mumps.nprow = 1
  mumps.npcol = 1
  mumps.mblock = 100
  mumps.nblock = 100
  mumps.schur_lld = mumps.size_schur
  schur = Array{T}(undef, mumps.size_schur^2)
  mumps.schur = pointer(schur)
  set_icntl!(mumps, 19, 3)
  append!(mumps._gc_haven, [Ref(listvar_schur), Ref(schur)])
  return mumps
end

function LinearAlgebra.det(mumps::Mumps{T}) where {T}
  if !has_det(mumps)
    throw(MUMPSException("ICNTL[33]=0, determinant not computed"))
  else
    if T <: Complex
      d = complex(mumps.rinfog[12], mumps.rinfog[13]) * T(2)^mumps.infog[34]
    else
      d = mumps.rinfog[12] * T(2)^mumps.infog[34]
    end
  end
  return T(d)
end

####################################################################
### auxilliary functions
####################################################################
is_finalized(mumps::Mumps) = mumps._finalized

is_matrix_assembled(mumps::Mumps) = !(mumps.icntl[5] ∈ [1])

is_matrix_distributed(mumps::Mumps) = mumps.icntl[18] ∈ [1, 2, 3]

is_rhs_dense(mumps::Mumps) = mumps.icntl[20] ∉ [1, 2, 3]

is_sol_central(mumps::Mumps) = mumps.icntl[21] ∉ [1]

has_det(mumps::Mumps) = mumps.icntl[33] ∉ [0]

is_symmetric(mumps::Mumps) = mumps.sym ∈ [1, 2]

is_posdef(mumps::Mumps) = mumps.sym ∈ [1]

has_matrix(mumps::Mumps) = mumps.n > 0

has_rhs(mumps::Mumps) = mumps.nrhs * mumps.lrhs > 0 || mumps.nz_rhs > 0

has_schur(mumps::Mumps) = mumps.size_schur > 0

LinearAlgebra.issymmetric(mumps::Mumps) = is_symmetric(mumps)
LinearAlgebra.isposdef(mumps::Mumps) = is_posdef(mumps)
