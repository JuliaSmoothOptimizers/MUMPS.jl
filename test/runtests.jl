using Base.Test
using MUMPS
using MPI

function random_matrix(t, σ, nrow, ncol)
  # t = data type
  # σ = vector of desired singular values
  Σ = diagm(Array{t}(σ))
  A = rand(t, nrow, ncol)
  (U, _) = qr(A)
  A = rand(t, nrow, ncol)
  (V, _) = qr(A)
  return U * Σ * V'
end

include("get_div_grad.jl")
root = 0

# Initialize MPI.
MPI.Init()
comm = MPI.COMM_WORLD

include("mumps_test_float.jl")
include("mumps_test.jl")
include("mumps_test_float_complex.jl")
include("mumps_test_complex.jl")
include("mumps_test_mixed.jl")

MPI.Barrier(comm)
MPI.Finalize()
