module MUMPS

export default_icntl, default_cntl, Mumps, finalize, factorize, solve,
       mumps_unsymmetric, mumps_definite, mumps_symmetric

# libjmumps.dylib should be on your LD_LIBRARY_PATH.
mumps_lib = "libjmumps";
macro mumps_call(func, args...)
  quote
    ccall(($func, $mumps_lib), $(args...))
  end
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

# Symbols for symmetry
mumps_unsymmetric = 0;
mumps_definite    = 1;
mumps_symmetric   = 2;


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
  det     :: Int
  err     :: Int

  # Constructor.
  function Mumps(sym   :: Int=mumps_unsymmetric;
                 icntl :: Array{Int32,1}=default_icntl,
                 cntl  :: Array{Float64,1}=default_cntl)
    # sym = 0 (unsymmetric), 1 (symmetric definite), 2 (symmetric indefinite).
    # There is currently no facility to exploit hermicity of complex matrices.

    # Set default pivot threshold if required.
    if cntl[1] == -1
      cntl[1] = (sym == mumps_definite) ? 0.0 : 0.01
    end

    id = @mumps_call(:mumps_initialize, Ptr{Void},
                     (Int32, Ptr{Int32}, Ptr{Float64}), sym, icntl, cntl);

    if id == C_NULL
      msg = "Error allocating MUMPS structure"
      error(msg)
    end

    infog = zeros(Int32, 40);
    rinfog = zeros(Float64, 20);

    self = new(id, int32(sym), icntl, cntl, 0, Float64, infog, rinfog, 0, 0, 0);
    finalizer(self, finalize);  # Destructor.
    return self;
  end
end


# Convenience constructor.
# Note that every argument is a keyword.
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
  if mumps.icntl[33] == 1
    mumps.det = mumps.rinfog[12] * 2^(mumps.infog[34]);
  end
  mumps.err = mumps.infog[1];
  return;
end


factorize(mumps :: Mumps, A :: Array{Float64}) = factorize(mumps, sparse(A));


function solve(mumps :: Mumps, rhs :: Array{Float64}; transposed :: Bool=false)
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
              (Ptr{Void}, Int32, Ptr{Float64},              Int32),
              mumps.__id,  nrhs,            x, transposed ? 1 : 0);

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
function solve(A :: SparseMatrixCSC, rhs :: Array{Float64};
               sym :: Int=mumps_unsymmetric)

  mumps = Mumps(sym);
  x = solve(mumps, A, rhs);
  return x;
end

end  # Module MUMPS
