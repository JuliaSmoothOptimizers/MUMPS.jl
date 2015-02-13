using Base.Test
using MUMPS
using MPI

include("get_div_grad.jl");
root = 0;

# Initialize MPI.
MPI.Init()
comm = MPI.COMM_WORLD

include("mumps_test.jl")
include("mumps_test_complex.jl")

MPI.Barrier(comm)
MPI.Finalize()

