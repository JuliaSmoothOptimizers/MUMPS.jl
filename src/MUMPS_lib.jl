"""Register the matrix `A` with the `Mumps` object `mumps`.
This function makes it possible to define the matrix on the host
only. If the matrix is defined on all nodes, there is no need to
use this function."""
function associate_matrix!{Tv <: MUMPSValueDataType, Ti <: MUMPSIntDataType}(mumps :: Mumps{Tv}, A :: SparseMatrixCSC{Tv,Ti})

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
associate_matrix!{Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti}) = associate_matrix!(mumps, convert(SparseMatrixCSC{Tm,Int64}, A));

# associate_matrix for dense matrices.
associate_matrix!{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, A :: Array{Tv,2}) = associate_matrix!(mumps, convert(SparseMatrixCSC{Tm,Int64}, sparse(A)));


# import Base.LinAlg.factorize

"""Factorize the matrix registered with the `Mumps` instance.
The matrix must have been previously registered with `associate_matrix()`.
After the factorization, the determinant, if requested, is stored in
`mumps.det`. The MUMPS error code is stored in `mumps.err`. """
function factorize!{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv})

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
on all nodes, there is no need to use this function.

This variant makes a copy of the right-hand side to avoid overwriting it."""
function associate_rhs{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}, rhs :: Array{Tv})
  associate_rhs!(mumps, copy(rhs))
end


"""Register the right-hand side(s) `rhs` with the `Mumps`
object `mumps`. This function makes it possible to define the right-
-hand side(s) on the host only. If the right-hand side(s) are defined
on all nodes, there is no need to use this function."""
function associate_rhs!{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}, rhs :: Array{Tv})

  n = size(rhs, 1);
  n == mumps.n || throw(MUMPSException("rhs has incompatible dimension"))

  nrhs = size(rhs, 2);
  x = rhs[:];  # Make a copy; will be overwritten with solution.

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
associate_rhs!{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, rhs :: Array{Tv}) = associate_rhs!(mumps, convert(Array{Tm}, rhs));


"""Solve the system registered with the `Mumps` object `mumps`.
The matrix and right-hand side(s) must have been previously registered
with `associate_matrix()` and `associate_rhs()`. The optional keyword
argument `transposed` indicates whether the user wants to solve the
forward or transposed system. The solution is stored internally and must
be retrieved with `get_solution()`."""
function solve!{Tv <: MUMPSValueDataType}(mumps :: Mumps{Tv}; transposed :: Bool=false)

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

  return reshape(x, int(mumps.n), nrhs);
end


# Convenience functions.

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
function factorize!{Tm <: MUMPSValueDataType, Tv <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti})
  mumps = associate_matrix!(mumps, A);  # A will be converted by associate_matrix.
  return factorize!(mumps);
end

"""Combined associate_matrix / factorize.
Presume that `A` is available on all nodes."""
factorize!{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, A :: Array{Tv}) = factorize!(mumps, convert(SparseMatrixCSC{Tm,Int64}, sparse(A)));


"""Combined associate_rhs / solve.
Presume that `rhs` is available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve{Tm <: MUMPSValueDataType, Tv <: Number}(mumps :: Mumps{Tm}, rhs :: Array{Tv}; transposed :: Bool=false)
  mumps = associate_rhs!(mumps, rhs);  # rhs will be converted by associate_rhs.
  mumps = solve!(mumps, transposed=transposed);
  return get_solution(mumps);
end


"""Combined analyze / factorize / solve.
Presume that `A` and `rhs` are available on all nodes.
The optional keyword argument `transposed` indicates whether
the user wants to solve the forward or transposed system.
The solution is retrieved and returned."""
function solve{Tm <: MUMPSValueDataType, Tv <: Number, Tr <: Number, Ti <: Integer}(mumps :: Mumps{Tm}, A :: SparseMatrixCSC{Tv,Ti}, rhs :: Array{Tr}; transposed :: Bool=false)

  mumps = factorize!(mumps, A);
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
