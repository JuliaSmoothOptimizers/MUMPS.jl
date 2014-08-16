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

Whether it uses the Accelerate framework or OpenBLAS, MUMPS links with 32 bit
versions of the BLAS and LAPACK libraries.

**You must make sure that Julia also uses 32 bit versions of the BLAS and
LAPACK**. In particular, if you use Homebrew to build Julia, **do not use** the
`--64` command-line argument to `brew install julia`.

At the Julia prompt, type

````JULIA
julia> Pkg.clone("https://github.com/dpo/MUMPS.jl.git")
````

Once MUMPS is cloned, building the interface is a matter of
````
cd ~/.julia/v0.x/MUMPS/src  # replace version number
make
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH  # Put this in your ~/.bashrc
````

## How to Use

The main data type holding information on a factorization is `Mumps`. Remember
to initialize MPI before attempting to create a `Mumps` object. A simple
session is as follows:

````JULIA
julia> using MUMPS
julia> ierr = mumps_initialize_mpi();
julia> A = sprand(10, 10, .2) + speye(10);
julia> rhs = rand(10);
julia> x = solve(A, rhs);
julia> norm(x - A \ rhs) / norm(x)
2.640677159735313e-16
julia> ierr = mumps_finalize_mpi();  # if you're finished
````

It is possible to separate the initialization, the analysis/factorization,
and the solve phases. It is also possible to access the information reported by
MUMPS after the factorization and solve phases, and to modify this information
(e.g., to perform iterative refinement). See the test examples.

## To Do (Pull Requests Welcome!)

* [ ] Support complex arithmetic
* [ ] Support distributed matrices / vectors

[![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)](http://www.gnu.org/licenses/gpl.html "GPLv3")
