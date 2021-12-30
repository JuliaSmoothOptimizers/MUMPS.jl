for (fname, elty) in (
  (:mumps_associate_matrix_float, Float32),
  (:mumps_associate_matrix_double, Float64),
  (:mumps_associate_matrix_float_complex, ComplexF32),
  (:mumps_associate_matrix_double_complex, ComplexF64),
)
  @eval begin
    """Register the matrix `A` with the `Mumps` object `mumps`.
    This function makes it possible to define the matrix on the host
    only. If the matrix is defined on all nodes, there is no need to
    use this function."""
    function associate_matrix!(
      mumps::Mumps{$elty},
      A::SparseMatrixCSC{$elty, Ti},
    ) where {Ti <: MUMPSIntDataType}
      n = size(A, 1)
      size(A, 2) == n || throw(MUMPSException("Input matrix must be square"))

      # Symmetric factorization only accesses the lower triangle.
      B = mumps.__sym > 0 ? tril(A) : A

      # Obtain B in coordinate format.
      nz = nnz(B)
      vals = convert(Array{$elty, 1}, B.nzval)    # Necessary?
      irow = convert(Array{Int32, 1}, B.rowval)   # Necessary?
      jcol = zeros(Int32, nz, 1)
      for i = 1:n
        jcol[B.colptr[i]:(B.colptr[i + 1] - 1)] .= i
      end

      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call(
        $(string(fname)),
        Nothing,
        (Ptr{Nothing}, Int32, Int32, Ptr{$elty}, Ptr{Int32}, Ptr{Int32}),
        id,
        n,
        nz,
        vals,
        irow,
        jcol
      )
      mumps.n = n
      mumps.nnz = mumps.infog[29]
      return mumps
    end

    function associate_matrix!(
      mumps::Mumps{$elty},
      n::Ti,
      irow::Vector{Ti},
      jcol::Vector{Ti},
      vals::Vector{$elty},
    ) where {Ti <: MUMPSIntDataType}
      nz = length(vals)
      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call(
        $(string(fname)),
        Nothing,
        (Ptr{Nothing}, Int32, Int32, Ptr{$elty}, Ptr{Int32}, Ptr{Int32}),
        id,
        n,
        nz,
        vals,
        irow,
        jcol
      )
      mumps.n = n
      mumps.nnz = mumps.infog[29]
      return mumps
    end
  end
end

# Associate a generally-typed matrix with a Mumps type. Attempt conversion.
# An InexactError should be raised if, e.g., mumps is Float64 and A is ComplexF64.
associate_matrix!(
  mumps::Mumps{Tm},
  A::SparseMatrixCSC{Tv, Ti},
) where {Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer} =
  associate_matrix!(mumps, convert(SparseMatrixCSC{Tm, Int32}, A))

# associate_matrix for dense matrices.
associate_matrix!(
  mumps::Mumps{Tm},
  A::Array{Tv, 2},
) where {Tm <: MUMPSValueDataType, Tv <: Number} =
  associate_matrix!(mumps, convert(SparseMatrixCSC{Tm, Int32}, sparse(A)))

for (fname, infoname, elty, infoty) in (
  (:mumps_factorize_float, :mumps_get_info_float, Float32, Float32),
  (:mumps_factorize_double, :mumps_get_info_double, Float64, Float64),
  (:mumps_factorize_float_complex, :mumps_get_info_float_complex, ComplexF32, Float32),
  (:mumps_factorize_double_complex, :mumps_get_info_double_complex, ComplexF64, Float64),
)
  @eval begin
    """Factorize the matrix registered with the `Mumps` instance.
    The matrix must have been previously registered with `associate_matrix()`.
    After the factorization, the determinant, if requested, is stored in
    `mumps.det`. The MUMPS error code is stored in `mumps.err`. """
    function factorize!(mumps::Mumps{$elty})
      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call($(string(fname)), Nothing, (Ptr{Nothing},), id)
      @mumps_call(
        $(string(infoname)),
        Nothing,
        (Ptr{Nothing}, Ptr{Int32}, Ptr{$infoty}),
        id,
        mumps.infog,
        mumps.rinfog
      )

      if mumps.icntl[33] == 1
        mumps.det = mumps.rinfog[12] * 2.0^(mumps.infog[34])
      end
      mumps.err = mumps.infog[1]
      return mumps
    end
  end
end

for (fname, elty) in (
  (:mumps_associate_rhs_float, Float32),
  (:mumps_associate_rhs_double, Float64),
  (:mumps_associate_rhs_float_complex, ComplexF32),
  (:mumps_associate_rhs_double_complex, ComplexF64),
)
  @eval begin
    """Register the right-hand side(s) `rhs` with the `Mumps`
    object `mumps`. This function makes it possible to define the right-
    -hand side(s) on the host only. If the right-hand side(s) are defined
    on all nodes, there is no need to use this function."""
    function associate_rhs!(mumps::Mumps{$elty}, rhs::Array{$elty})
      n = size(rhs, 1)
      n == mumps.n || throw(MUMPSException("rhs has incompatible dimension"))

      nrhs = size(rhs, 2)
      x = rhs[:]  # Make a copy; will be overwritten with solution.

      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call($(string(fname)), Nothing, (Ptr{Nothing}, Int32, Ptr{$elty}), id, nrhs, x)
      return mumps
    end
  end
end

# Associate a generally-typed rhs with a Mumps type. Attempt conversion.
# An InexactError should be raised if, e.g., mumps is Float64 and rhs is ComplexF64.
associate_rhs!(mumps::Mumps{Tm}, rhs::Array{Tv}) where {Tm <: MUMPSValueDataType, Tv <: Number} =
  associate_rhs!(mumps, convert(Array{Tm}, rhs))

for (fname, infoname, elty, infoty) in (
  (:mumps_solve_float, :mumps_get_info_float, Float32, Float32),
  (:mumps_solve_double, :mumps_get_info_double, Float64, Float64),
  (:mumps_solve_float_complex, :mumps_get_info_float_complex, ComplexF32, Float32),
  (:mumps_solve_double_complex, :mumps_get_info_double_complex, ComplexF64, Float64),
)
  @eval begin
    """Solve the system registered with the `Mumps` object `mumps`.
    The matrix and right-hand side(s) must have been previously registered
    with `associate_matrix()` and `associate_rhs()`. The optional keyword
    argument `transposed` indicates whether the user wants to solve the
    forward or transposed system. The solution is stored internally and must
    be retrieved with `get_solution()`."""
    function solve!(mumps::Mumps{$elty}; transposed::Bool = false)
      id = reinterpret(Ptr{Nothing}, mumps.__id)
      @mumps_call($(string(fname)), Nothing, (Ptr{Nothing}, Int32), id, transposed ? 1 : 0)
      @mumps_call(
        $(string(infoname)),
        Nothing,
        (Ptr{Nothing}, Ptr{Int32}, Ptr{$infoty}),
        id,
        mumps.infog,
        mumps.rinfog
      )

      mumps.err = mumps.infog[1]
      return mumps
    end
  end
end

for (fname, solname, elty) in (
  (:mumps_get_nrhs_float, :mumps_get_solution_float, Float32),
  (:mumps_get_nrhs_double, :mumps_get_solution_double, Float64),
  (:mumps_get_nrhs_float_complex, :mumps_get_solution_float_complex, ComplexF32),
  (:mumps_get_nrhs_double_complex, :mumps_get_solution_double_complex, ComplexF64),
)
  @eval begin
    """Retrieve the solution of the system solved by `solve()`. This
    function makes it possible to ask MUMPS to assemble the final solution
    on the host only, and to retrieve it there."""
    function get_solution(mumps::Mumps{$elty})
      id = reinterpret(Ptr{Nothing}, mumps.__id)
      nrhs = Int(@mumps_call($(string(fname)), Int32, (Ptr{Nothing},), id))
      x = zeros($elty, mumps.n * nrhs)
      @mumps_call($(string(solname)), Nothing, (Ptr{Nothing}, Ptr{$elty}), id, x)

      return reshape(x, Int(mumps.n), nrhs)
    end
  end
end

# Convenience functions.

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
function factorize!(
  mumps::Mumps{Tm},
  A::SparseMatrixCSC{Tv, Ti},
) where {Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer}
  mumps = associate_matrix!(mumps, A)  # A will be converted by associate_matrix.
  return factorize!(mumps)
end

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
factorize!(mumps::Mumps{Tm}, A::Array{Tv}) where {Tm <: MUMPSValueDataType, Tv <: Number} =
  factorize!(mumps, convert(SparseMatrixCSC{Tm, Int32}, sparse(A)))

"""Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve(
  mumps::Mumps{Tm},
  rhs::Array{Tv};
  transposed::Bool = false,
) where {Tm <: MUMPSValueDataType, Tv <: Number}
  mumps = associate_rhs!(mumps, rhs)  # rhs will be converted by associate_rhs.
  mumps = solve!(mumps, transposed = transposed)
  return get_solution(mumps)
end

"""Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve(
  mumps::Mumps{Tm},
  A::SparseMatrixCSC{Tv, Ti},
  rhs::Array{Tr};
  transposed::Bool = false,
) where {Tm <: MUMPSValueDataType, Tv <: Number, Tr <: Number, Ti <: Integer}
  mumps = factorize!(mumps, A)
  return solve(mumps, rhs, transposed = transposed)
end

solve(
  mumps::Mumps{Tm},
  A::Array{Tv, 2},
  rhs::Array{Tr},
) where {Tm <: MUMPSValueDataType, Tv <: Number, Tr <: Number} =
  solve(mumps, convert(SparseMatrixCSC{Tv, Int32}, sparse(A)), rhs)

"""Combined initialize / analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `sym` indicates the symmetry of `A`.
The solution is retrieved and returned."""
function solve(
  A::SparseMatrixCSC{Tv, Ti},
  rhs::Array{Tr};
  sym::Int = mumps_unsymmetric,
) where {Tv <: Number, Tr <: Number, Ti <: MUMPSIntDataType}
  Tm = (Tv <: Complex || Tr <: Complex) ? ComplexF64 : Float64  # Could be smarter.
  mumps = Mumps{Tm}(sym, default_icntl, default_cntl64)
  x = solve(mumps, A, rhs)
  finalize(mumps)
  return x
end

solve(
  A::Array{Tv, 2},
  rhs::Array{Tr};
  sym::Int = mumps_unsymmetric,
) where {Tv <: Number, Tr <: Number} =
  solve(convert(SparseMatrixCSC{Tv, Int32}, sparse(A)), rhs, sym = sym)
