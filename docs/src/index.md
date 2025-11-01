# [MUMPS.jl documentation](@id Home)

MUMPS is a library for the solution of large linear systems using a
factorization. Structure can be exploited, such as symmetry, or symmetry and
definiteness. The factorization and solve phases can be performed in parallel
via MPI by way of [MPI.jl](https://github.com/JuliaParallel/MPI.jl).

## How to Install

```julia
julia> ]
pkg> add MUMPS
pkg> test MUMPS
```

## How to Use

The main data type holding information on a factorization is `Mumps`. Remember
to initialize MPI before attempting to create a `Mumps` object. A simple
session is as follows:

```julia
using MUMPS, MPI, SparseArrays, LinearAlgebra
MPI.Init()
A = sprand(10, 10, 0.2) + I
rhs = rand(10)
x = solve(A, rhs)
norm(x - A \ rhs) / norm(x)
MPI.Finalize()
```

It is possible to separate the initialization, the analysis/factorization,
and the solve phases. It is also possible to access the information reported by
MUMPS after the factorization and solve phases, and to modify this information
(e.g., to perform iterative refinement).

When creating an instance of a `Mumps` object explicitly, it is important to
specify in advance what arithmetic should be used. Single and double precision
real (`Float32` and `Float64`) and complex (`ComplexF32` and `ComplexF64`)
arithmetics are supported.

For instance,

```julia
MPI.Init()
mumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64)
A = sparse(rand(4,4))
rhs = rand(4)
associate_matrix!(mumps, A)
factorize!(mumps)
associate_rhs!(mumps, rhs)
solve!(mumps)
x = get_solution(mumps)
finalize(mumps)
MPI.Finalize()
```

Once the arithmetic of the `Mumps` instance has been specified, it cannot be
changed. The module is flexible in that various data types may be used to
define the matrix to be factorized and the right-hand side, and appropriate
conversions will take place. Dense matrices may be used, and they will be
converted to sparse format.

For intance,

```julia
mumps = Mumps{ComplexF64}(mumps_unsymmetric, default_icntl, default_cntl64)
A = rand(Int16, 4, 4)
rhs = rand(Float32, 4)
associate_matrix!(mumps, A)  # A is converted to a sparse ComplexF64 matrix
associate_rhs!(mumps, rhs)   # rhs is converted to a Complex64 vector
```

See [test](https://github.com/JuliaSmoothOptimizers/MUMPS.jl/tree/main/test) for more examples.

## Constants and Methods Exposed

### Constants

The following convenience constants may be used when initializing a `Mumps`
object:

Constant            | Meaning
--------------------|--------
`mumps_unsymmetric` | matrix is general unsymmetric (or symmetry is unknown)
`mumps_definite`    | matrix is symmetric and (positive or negative) definite
`mumps_symmetric`   | matrix is symmetric but indefinite (or definiteness is unknown)
`default_icntl`     | array of default integer parameters
`default_cntl32`    | array of default real parameters in single precision
`default_cntl64`    | array of default real parameters in double precision

See Section 5 of the [MUMPS User's Manual](https://mumps-solver.org/index.php?page=doc) for a description of the integer and real control arrays.

### Methods

A `Mumps` object is created using the default constructor, which must be
supplied with:

* the data type for the arithmetic to be used, as a type parameter, i.e.,
  `Mumps{Float64}(...)` or `Mumps{ComplexF64}(...)`
* `sym`: one of the constants `mumps_unsymmetric`, `mumps_definite` or
  `mumps_symmetric`. Note that there is no support for Hermitian complex
  matrices in MUMPS. Therefore, we recommend to always use `mumps_unsymmetric`
  for complex data.
* `icntl`: an integer parameters array (see the MUMPS Users's Manual)
* `cntl`: a real parameters array (see the MUMPS Users's Manual)

The convenience function `get_icntl()` returns an array of integer parameters
corresponding to certain commonly-used options. Its arguments are all optional:

* `det`: a boolean indicating whether the determinant should be computed
  (default: `false`)
* `verbose`: a boolean (default: `false`)
* `ooc`: a boolean indicating whether factors should be stored out of core
  (default: `false`)
* `itref`: the number of iterative refinement steps (default: 0).

A `Mumps` object is destroyed by calling the `finalize()` method. Because
`finalize` still issues MPI commands, it is important to call `finalize()`
before calling `MPI.Finalize()`.

Method                    | Description
--------------------------|------------
`finalize`                | Finalize a `Mumps` object. Must be done before calling `MPI.Finalize()`
`associate_matrix!`       | Register a matrix with the `Mumps` object. This function makes it possible to define the data on the host only.
`factorize!`              | Factorize the matrix registered with the `Mumps` object.
`associate_rhs!`          | Register right-hand sides with the `Mumps` object. This function makes it possible to define the data on the host only.
`solve!`                  | Solve the linear system for the given right-hand side.
`get_solution`            | Retrieve the solution from the `Mumps` object. This function makes it possible for the solution to be assembled on the host only.

## Parallel Execution

MPI is controled by way of [MPI.jl](https://github.com/JuliaParallel/MPI.jl).
Look for the lines that say `NUMBER OF WORKING PROCESSES` in the output of

```shell
mpirun -np 4 julia examples/mumps_mpi.jl
```

## Distributed arrays (DArray)

`DistributedArrays` has been added for workflows that use distributed matrices and vectors. It allows tests to validate behavior when data is stored in `DArray` form and interoperates with MPI-backed workflows.

```julia
using Distributed, DistributedArrays, MPI, MUMPS, SparseArrays, LinearAlgebra

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

sqrtN = 40
N = sqrtN^2

function local_rows(I)
  # DistributedArrays may pass a 1-tuple of a UnitRange to the init function on workers; extract the actual index range if necessary.
  idx = isa(I, Tuple) ? I[1] : I
  nn = length(idx)
  # Return a vector of sparse row vectors, one per local row. This makes the DArray hold per-row SparseVector objects so the root can reassemble the full sparse matrix from triplets.
  rows_vec = Vector{SparseVector{Float64}}(undef, nn)
  for (local_i, global_i) in enumerate(idx)
    i = div(global_i-1, sqrtN) + 1
    j = (global_i-1) % sqrtN + 1
    cols = Int[]; vals = Float64[]
    push!(cols, global_i); push!(vals, 4.0)
    if i > 1  push!(cols, global_i - sqrtN); push!(vals, -1.0) end
    if i < sqrtN push!(cols, global_i + sqrtN); push!(vals, -1.0) end
    if j > 1  push!(cols, global_i - 1); push!(vals, -1.0) end
    if j < sqrtN push!(cols, global_i + 1); push!(vals, -1.0) end
    rows_vec[local_i] = sparsevec(cols, vals, N)
  end
  return rows_vec
end

# Create a DArray of local sparse row-blocks
dblocks = DArray(I -> local_rows(I), (N,))

if rank == 0
  # Gather local blocks and assemble global sparse matrix on the root.
  # Use vertical concatenation of sparse blocks to avoid slice-assignment edge
  # cases when assigning sparse blocks into a sparse matrix.
  # Collect per-row sparse vectors from the DArray and build triplets
  blocks = collect(dblocks)   # blocks is a Vector of SparseVector objects
  rows_idx = Int[]; cols_idx = Int[]; vals = Float64[]
  global_row = 1
  for sv in blocks
    if nnz(sv) > 0
      inds, valsv = findnz(sv)
      append!(rows_idx, fill(global_row, length(inds)))
      append!(cols_idx, inds)
      append!(vals, valsv)
    end
    global_row += 1
  end
  A = sparse(rows_idx, cols_idx, vals, N, N)

  b = ones(Float64, N)

  # Enable ScaLAPACK usage on the root frontal matrix
  icntl = default_icntl[:]
  icntl[13] = 0
  m = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  associate_matrix!(m, A)
  factorize!(m)
  associate_rhs!(m, b)
  solve!(m)
  x = get_solution(m)
  println("Residual norm on root: ", norm(A * x - b))
  finalize(m)
else
  # non-root ranks finished after local assembly
end

MPI.Barrier(comm)
MPI.Finalize()
```

### ScaLAPACK Support

MUMPS is compiled with **ScaLAPACK** and **PARMETIS** support, which provides improved performance for parallel factorization, particularly for operations involving the Schur complement on the root node when using MPI parallelism.

The use of ScaLAPACK is controlled by the `ICNTL[13]` parameter:

* **`ICNTL[13] = 0`** (default): Enables parallel factorization using ScaLAPACK on the root frontal matrix. This is the recommended setting for most parallel computations.
* **`ICNTL[13] > 0`**: Forces sequential factorization on the root frontal matrix unless the number of working processes exceeds the specified value. Use this if you want to disable ScaLAPACK or control when it is used based on the number of workers.

#### Example: Controlling ScaLAPACK usage

```julia
using MUMPS, MPI, SparseArrays

MPI.Init()

icntl = default_icntl[:]
icntl[13] = 0

mumps = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
A = sprand(100, 100, 0.1)
factorize!(mumps, A)

finalize(mumps)
MPI.Finalize()
```

For more details on control parameters, see Section 5 of the [MUMPS User's Manual](https://mumps-solver.org/index.php?page=doc).

## Custom Installation

**Note: MUMPS is already precompiled with Yggdrasil for all platforms except Windows.**

To use your custom MUMPS, set the environment variable `JULIA_MUMPS_LIBRARY_PATH`
to point to the shared library before `using MUMPS`.
Note that the same version of MUMPS as used by the `MUMPS_jll` artifact is needed.

For example, macOS users may install precompiled MUMPS binaries from the Homebrew tap `dpo/mumps-jl` as follows:

```bash
brew tap dpo/mumps-jl
brew install mpich-mumps
export JULIA_MUMPS_LIBRARY_PATH=$(brew --prefix)/opt/mpich-mumps/lib
```

Apple Silicon users should remember to use `arch x86_64 brew` to refer to Intel binaries run through Rosetta, as we do not (yet) ship Silicon binaries of MUMPS via Homebrew.

The `JULIA_MUMPS_LIBRARY_PATH` environment variable may be set permanently in the shell's startup file, or in `$HOME/.julia/config/startup.jl`.
