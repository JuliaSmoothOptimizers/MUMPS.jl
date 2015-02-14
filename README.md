# A [Julia](http://julialang.org) Interface to [MUMPS](http://mumps.enseeiht.fr)

OSX and Linux: [![Build Status](https://travis-ci.org/dpo/MUMPS.jl.svg?branch=master)](https://travis-ci.org/dpo/MUMPS.jl)
[![Coverage Status](https://coveralls.io/repos/dpo/MUMPS.jl/badge.svg)](https://coveralls.io/r/dpo/MUMPS.jl)

MUMPS is a library for the solution of large linear systems using a
factorization. Structure can be exploited, such as symmetry, or symmetry and
definiteness. The factorization and solve phases can be performed in parallel.

## Documentation

See http://dpo.github.io/MUMPS.jl/MUMPS.html for preliminary documentation.

## How to Install

On OSX we recommend using [Homebrew](https://brew.sh). The procedure below
should also work on Linux if you use [Linuxbrew](https://github.com/Homebrew/linuxbrew).

````
brew tap homebrew/science  # if not already done
brew install mumps [--with-scotch5] [--with-openblas]  # use scotch/openblas at your option
````

At the Julia prompt, type

````JULIA
julia> Pkg.clone("https://github.com/dpo/MUMPS.jl.git")
julia> Pkg.build("MUMPS")
````

In order for Julia to find the MUMPS interface library, its location must
appear on your `LD_LIBRARY_PATH`:
````
export LD_LIBRARY_PATH=$(julia -E 'Pkg.dir()' | sed -e 's/"//g')/MUMPS/src:$LD_LIBRARY_PATH
````

Place the above in your `~/.bashrc` to make it permanent.

## How to Use

The main data type holding information on a factorization is `Mumps`. Remember
to initialize MPI before attempting to create a `Mumps` object. A simple
session is as follows:

````JULIA
julia> using MUMPS
julia> using MPI
julia> MPI.Init()
julia> A = sprand(10, 10, .2) + speye(10); rhs = rand(10);
julia> x = solve(A, rhs);  # Mumps object is created and destroyed
julia> norm(x - A \ rhs) / norm(x)
2.640677159735313e-16
julia> MPI.Finalize();     # if you're finished
````

It is possible to separate the initialization, the analysis/factorization,
and the solve phases. It is also possible to access the information reported by
MUMPS after the factorization and solve phases, and to modify this information
(e.g., to perform iterative refinement).

When creating an instance of a `Mumps` object explicitly, it is important to
specify in advance what arithmetic should be used. Single and double precision
real (`Float32` and `Float64`) and complex (`Complex64` and `Complex128`)
arithmetics are supported.

For instance,

````JULIA
julia> MPI.Init();
julia> mumps = Mumps{Float64}(mumps_unsymmetric);  # Real, general unsymmetric
julia> A = sparse(rand(4,4)); rhs = rand(4);       # Happens on all cores
julia> associate_matrix(mumps, A);
julia> factorize(mumps);
julia> associate_matrix(mumps, rhs);
julia> solve(mumps);
julia> x = get_solution(mumps);
julia> finalize(mumps);
julia> MPI.Finalize();
````

Once the arithmetic of the `Mumps` instance has been specified, it cannot be
changed. The module is flexible in that various data types may be used to
define the matrix to be factorized and the right-hand side, and appropriate
conversions will take place. Dense matrices may be used, and they will be
converted to sparse format.

For intance,

```JULIA
julia> mumps = Mumps{Complex128}(mumps_unsymmetric);
julia> A = rand(Int16, 4, 4); rhs = rand(Float32, 4);
julia> associate_matrix(mumps, A);  # A is converted to a sparse Float64 matrix
julia> associate_rhs(mumps, rhs);   # rhs is converted to a Float64 array
```

See [test](https://github.com/dpo/MUMPS.jl/tree/master/test) for more examples.

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

See Sections 5.1 and 5.2 of the [MUMPS User's Manual](http://mumps.enseeiht.fr/doc/userguide_4.10.0.pdf) for a description of the integer and real control arrays.

### Methods

A `Mumps` object is created using the default constructor, which must be
supplied with:

* the data type for the arithmetic to be used, as a type parameter, i.e.,
  `Mumps{Float64}(...)` or `Mumps{Complex128}(...)`
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
`associate_matrix` | Register a matrix with the `Mumps` object. This function makes it possible to define the data on the host only.
`factorize`        | Factorize the matrix registered with the `Mumps` object.
`associate_rhs`    | Register right-hand sides with the `Mumps` object. This function makes it possible to define the data on the host only.
`solve`            | Solve the linear system for the given right-hand side.
`get_solution`     | Retrieve the solution from the `Mumps` object. This function makes it possible for the solution to be assembled on the host only.

## Parallel Execution

MPI is controled by way of [MPI.jl](https://github.com/lcw/MPI.jl).
Look for the lines that say `NUMBER OF WORKING PROCESSES` in the output of

````
mpirun -np 4 julia examples/mumps_mpi.jl
````

## To Do (Pull Requests Welcome!)

* [X] Support double precision complex arithmetic (in [99c23fe](https://github.com/dpo/MUMPS.jl/commit/99c23fe87e7c985fe3062d78ab7664b82a6b8dba))
* [X] Support single precision real and complex arithmetic (in [654814a](https://github.com/dpo/MUMPS.jl/commit/654814a5e5800260011d2f26f7fb6de179609cfa))
* [ ] Support distributed matrices / vectors
* [ ] User-selected permutation
* [X] Out-of-core option (in [73e829b](https://github.com/dpo/MUMPS.jl/commit/73e829b52fe3d20c70c2733607ba9820cda03ed6#diff-d41d8cd98f00b204e9800998ecf8427e))
* [X] Determinant (in [73e829b](https://github.com/dpo/MUMPS.jl/commit/73e829b52fe3d20c70c2733607ba9820cda03ed6#diff-d41d8cd98f00b204e9800998ecf8427e))
* [ ] Compute entries of the inverse
* [X] Control iterative refinement (in [73e829b](https://github.com/dpo/MUMPS.jl/commit/73e829b52fe3d20c70c2733607ba9820cda03ed6#diff-d41d8cd98f00b204e9800998ecf8427e))
* [ ] Obtain a Schur complement
* [ ] Solve with sparse right-hand sides

This content is released under the [MIT](http://opensource.org/licenses/MIT) License.
<a rel="license" href="http://opensource.org/licenses/MIT">
<img alt="MIT license" height="40" src="http://upload.wikimedia.org/wikipedia/commons/c/c3/License_icon-mit.svg" /></a>
