#!/bin/bash
set -ev

export JULIA_PKG_DIR=$(julia -E 'Pkg.dir()' | sed -e 's/"//g')
export LD_LIBRARY_PATH=$JULIA_PKG_DIR/MUMPS/src:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
ls -l $JULIA_PKG_DIR/
ls -l $JULIA_PKG_DIR/MUMPS/src
julia --check-bounds=yes -E 'Pkg.test("MUMPS")'
mpirun -np 2 julia $TRAVIS_BUILD_DIR/examples/mumps_mpi.jl
