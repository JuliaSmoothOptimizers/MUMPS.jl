# A [Julia](http://julialang.org) Interface to [MUMPS](http://mumps.enseeiht.fr)

MUMPS is a library for the solution of large linear systems using a
factorization. Structure can be exploited, such as symmetry, or symmetry and
definiteness. The factorization and solve phases can be performed in parallel.

## How to Install

We recommend using [Homebrew](https://brew.sh). Currently, only OSX is
supported. The procedure below may work if you use
[Linuxbrew](https://github.com/Homebrew/linuxbrew).

````
brew tap homebrew/science  # if not already done
brew install mumps --with-mpi --with-scotch5 [--with-openblas]  # use openblas at your option
````

At the Julia prompt, type

````JULIA
julia> Pkg.clone("https://github.com/dpo/MUMPS.jl.git")
julia> Pkg.build("MUMPS")
````

In order for Julia to find the MUMPS interface library, its location must
appear on your `LD_LIBRARY_PATH`:
````
export LD_LIBRARY_PATH=~/.julia/v0.x/MUMPS/src:$LD_LIBRARY_PATH
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
