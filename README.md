# A [Julia](http://julialang.org) Interface to [MUMPS](http://mumps.enseeiht.fr)

OSX: [![Build
Status](https://travis-ci.org/dpo/MUMPS.jl.svg?branch=master)](https://travis-ci.org/dpo/MUMPS.jl)
Linux: Coming up...

MUMPS is a library for the solution of large linear systems using a
factorization. Structure can be exploited, such as symmetry, or symmetry and
definiteness. The factorization and solve phases can be performed in parallel.

## Documentation

See http://dpo.github.io/MUMPS.jl/MUMPS.html for preliminary documentation.

## How to Install

We recommend using [Homebrew](https://brew.sh). Currently, only OSX is
supported. The procedure below may work if you use
[Linuxbrew](https://github.com/Homebrew/linuxbrew).

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
(e.g., to perform iterative refinement). For instance,

````JULIA
julia> MPI.Init();
julia> mumps = Mumps(mumps_unsymmetric);  # General unsymmetric
julia> A = sparse(rand(4,4)); rhs = rand(4);  # Happens on all cores
julia> associate_matrix(mumps, A);
julia> factorize(mumps);
julia> associate_matrix(mumps, rhs);
julia> solve(mumps);
julia> x = get_solution(mumps);
julia> finalize(mumps);
julia> MPI.Finalize();
````

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
`default_cntl`      | array of default real parameters


See Sections 5.1 and 5.2 of the [MUMPS User's Manual](http://mumps.enseeiht.fr/doc/userguide_4.10.0.pdf) for a description of the integer and real control arrays.

### Methods

A `Mumps` object can be created in two ways

1. The convenience constructor has no required argument, but may optionally
   be supplied with the following keyword arguments:

    * `sym`: one of the constants `mumps_unsymmetric`, `mumps_definite` or `mumps_symmetric`
    * `det`: a boolean indicating whether the determinant should be computed
    * `verbose`: a boolean
    * `ooc`: a boolean indicating whether factors should stored out of core
    * `itref`: the number of iterative refinement steps
    * `cntl`: a real parameters array (see the MUMPS Users's Manual)

2. The standard constructor must be supplied with:

    * `sym`: one of the constants `mumps_unsymmetric`, `mumps_definite` or `mumps_symmetric`
    * `icntl`: an integer parameters array (see the MUMPS Users's Manual)
    * `cntl`: a real parameters array (see the MUMPS Users's Manual)

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

* [ ] Support complex arithmetic
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
