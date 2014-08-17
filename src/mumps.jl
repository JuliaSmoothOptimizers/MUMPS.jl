module MUMPS

export default_icntl, default_cntl, Mumps, finalize, factorize, solve,
       mumps_initialize_mpi, mumps_finalize_mpi

# libjmumps.dylib should be on your LD_LIBRARY_PATH.
mumps_lib = "libjmumps.dylib";
macro mumps_call(func, args...)
  quote
    ccall(($func, $mumps_lib), $(args...))
  end
end


# MPI-related functions.
function mumps_initialize_mpi()
  ierr = @mumps_call(:mumps_initialize_mpi, Int32, ());
  return ierr
end


function mumps_finalize_mpi()
  ierr = @mumps_call(:mumps_finalize_mpi, Int32, ())
  return ierr
end


# See MUMPS User's Manual Section 5.1.
default_icntl = convert(Array{Int32,1},
                        [6, 0,  6, 2, 0, 7,  7, 77, 1, 0,
                         0, 0, 0, 20, 0, 0,  0,  0, 0, 0,
                         0, 0, 0,  0, 0, 0, -8,  0, 0, 0,
                         0, 0, 0,  0, 0, 0,  0,  0, 0, 0]);

# See MUMPS User's Manual Section 5.2.
# icnt[1] will be set to its default value if left at -1.
default_cntl = [-1, sqrt(eps(1.0)), 0.0, -1.0, 0.0];


type Mumps
  __id    :: Ptr{Void}         # Pointer to MUMPS struct. Do not touch.
  __sym   :: Int32             # Value of sym used by Mumps.
  icntl   :: Array{Int32,1}    # Integer control parameters.
  cntl    :: Array{Float64,1}  # Real control parameters.
  n       :: Int32             # Order of the matrix factorized.
  valtype :: Union(Type{Float64}, Type{Complex64})
  infog   :: Array{Int32,1}
  rinfog  :: Array{Float64,1}
  nnz     :: Int               # Number of nonzeros in factors.
  err     :: Int

  # Constructor.
  function Mumps(sym :: Int=0,
                 icntl :: Array{Int32,1}=default_icntl,
                 cntl :: Array{Float64,1}=default_cntl)
    # sym = 0 (unsymmetric), 1 (symmetric definite), 2 (symmetric indefinite).
    # There is currently no facility to exploit hermicity of complex matrices.

    # Set default pivot threshold if required.
    if cntl[1] == -1
      cntl[1] = (sym > 0) ? 0.0 : 0.01
    end

    id = @mumps_call(:mumps_initialize, Ptr{Void},
                     (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);

    if id == C_NULL
      msg = "Error allocating MUMPS structure"
      error(msg)
    end

    infog = zeros(Int32, 40);
    rinfog = zeros(Float64, 20);

    self = new(id, int32(sym), icntl, cntl, 0, Float64, infog, rinfog, 0, 0);
    finalizer(self, finalize);  # Destructor.
    return self;
  end
end


function finalize(mumps :: Mumps)
  # Terminate a Mumps instance.
  @mumps_call(:mumps_finalize, Void, (Ptr{Void},), mumps.__id);
  mumps.__id = C_NULL;
end


function factorize(mumps :: Mumps, A :: SparseMatrixCSC)

  n = size(A, 1);
  if size(A, 2) != n
    error("Input matrix must be square");
  end

  # Symmetric factorization only accesses the lower triangle.
  B = mumps.__sym > 0 ? tril(A) : A;

  # Obtain B in coordinate format.
  nz = nnz(B);
  valtype = isreal(B.nzval[1]) ? Float64 : Complex64;
  vals = convert(Array{valtype,1}, B.nzval);
  irow = convert(Array{Int32,1}, B.rowval);
  jcol = zeros(Int32, nz, 1);
  for i = 1 : n
    jcol[B.colptr[i] : B.colptr[i+1]-1] = i;
  end

  @mumps_call(:mumps_factorize, Void,
              (Ptr{Void}, Int32, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
              mumps.__id,     n,    nz,         vals,       irow,       jcol);

  @mumps_call(:mumps_get_info, Void,
              (Ptr{Void}, Ptr{Int32}, Ptr{Float64}),
              mumps.__id, mumps.infog, mumps.rinfog)

  mumps.n = n;
  mumps.valtype = valtype;
  mumps.nnz = mumps.infog[29];
  mumps.err = mumps.infog[1];
  return;
end


function solve(mumps :: Mumps, rhs :: Array{Float64})
  n = size(rhs, 1);
  if n != mumps.n
    error("rhs has incompatible dimension");
  end
  if typeof(rhs[1]) != mumps.valtype
    error("rhs has incompatible data type");
  end

  nrhs = size(rhs, 2);
  x = rhs[:];
  @mumps_call(:mumps_solve, Void,
              (Ptr{Void}, Int32, Ptr{Float64}),
              mumps.__id,  nrhs,            x);

  @mumps_call(:mumps_get_info, Void,
              (Ptr{Void}, Ptr{Int32}, Ptr{Float64}),
              mumps.__id, mumps.infog, mumps.rinfog)

  mumps.err = mumps.infog[1];
  return reshape(x, n, nrhs);
end


# Convenience functions.

## Combined analyze / factorize / solve.
function solve(mumps :: Mumps, A :: SparseMatrixCSC, rhs :: Array{Float64})
  factorize(mumps, A);
  x = solve(mumps, rhs);
  return x;
end


## Combined initialize / analyze / factorize / solve.
function solve(A :: SparseMatrixCSC, rhs :: Array{Float64}, sym :: Int=0)
  mumps = Mumps(sym);  # Assume non-symmetric.
  x = solve(mumps, A, rhs);
  return x;
end

end  # Module MUMPS
