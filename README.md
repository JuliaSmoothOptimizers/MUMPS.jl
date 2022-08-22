# A [Julia](http://julialang.org) Interface to [MUMPS](https://graal.ens-lyon.fr/MUMPS/index.php?page=home)

[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![build-gh][build-gh-img]][build-gh-url]
[![build-cirrus][build-cirrus-img]][build-cirrus-url]
[![codecov][codecov-img]][codecov-url]
[![doi][doi-img]][doi-url]

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/MUMPS.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/MUMPS.jl/dev
[build-gh-img]: https://github.com/JuliaSmoothOptimizers/MUMPS.jl/workflows/CI/badge.svg?branch=main
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/MUMPS.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/JuliaSmoothOptimizers/MUMPS.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/JuliaSmoothOptimizers/MUMPS.jl
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/MUMPS.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/JuliaSmoothOptimizers/MUMPS.jl
[doi-img]: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.3271646-blue.svg
[doi-url]: https://doi.org/10.5281/zenodo.3271646

MUMPS is a library for the solution of large square linear systems using a
factorization. Structure can be exploited, such as symmetry, or symmetry and
definiteness. The factorization and solve phases can be performed in parallel
via MPI by way of [MPI.jl](https://github.com/JuliaParallel/MPI.jl).

## How to Cite

If you use MUMPS.jl in your work, please cite using the format given in [`CITATION.bib`](https://github.com/JuliaSmoothOptimizers/MUMPS.jl/blob/main/CITATION.bib)

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
using MUMPS, MPI, SparseArrays
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

See Section 5 of the [MUMPS User's Manual](https://graal.ens-lyon.fr/MUMPS/index.php?page=doc) for a description of the integer and real control arrays.

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

Method             | Description
-------------------|------------
`finalize`         | Finalize a `Mumps` object. Must be done before calling `MPI.Finalize()`
`associate_matrix!`| Register a matrix with the `Mumps` object. This function makes it possible to define the data on the host only.
`factorize!`       | Factorize the matrix registered with the `Mumps` object.
`associate_rhs!`   | Register right-hand sides with the `Mumps` object. This function makes it possible to define the data on the host only.
`solve!`           | Solve the linear system for the given right-hand side.
`get_solution`     | Retrieve the solution from the `Mumps` object. This function makes it possible for the solution to be assembled on the host only.

## Parallel Execution

MPI is controled by way of [MPI.jl](https://github.com/JuliaParallel/MPI.jl).
Look for the lines that say `NUMBER OF WORKING PROCESSES` in the output of

```shell
mpirun -np 4 julia examples/mumps_mpi.jl
```

## Custom Installation

**Note: MUMPS is already precompiled with Yggdrasil for all platforms except Windows.**

To use your custom MUMPS, set the environment variable `JULIA_MUMPS_LIBRARY_PATH`
to point to the shared library before `using MUMPS`. Note that **MUMPS** version 5.5.1 is needed.

For example:
```julia
ENV["JULIA_MUMPS_LIBRARY_PATH"] = "~/Applications/MUMPS_5.5.1/lib"
using MUMPS
```

Alternatively, you can set these permanently through your operating system.
