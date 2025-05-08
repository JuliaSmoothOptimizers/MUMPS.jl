# A [Julia](http://julialang.org) Interface to [MUMPS](https://mumps-solver.org/index.php?page=home)

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

MUMPS is a libray for the solution of sparse linear systems on multicore computers.
MUMPS implements methods based on the LDL or LU factorization of the input matrix and is suited to solving square symmetric or unsymmetric linear systems.
MUMPS supports real and complex, single and double precision arithmetic.

## How to Cite

If you use MUMPS.jl in your work, please cite using the format given in [`CITATION.cff`](https://github.com/JuliaSmoothOptimizers/MUMPS.jl/blob/main/CITATION.cff)
