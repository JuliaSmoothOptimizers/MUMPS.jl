using Test
using Random
using LinearAlgebra
using MPI
using MUMPS
using SparseArrays

Random.seed!(666)  # Random tests are diabolical

@info("MUMPS_INSTALLATION: $(MUMPS.MUMPS_INSTALLATION)")

function random_matrix(t, σ, nrow, ncol)
  # t = data type
  # σ = vector of desired singular values
  Σ = diagm(0 => Array{t}(σ))
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

@testset "float: " begin
  include("mumps_test_float.jl")
end
@testset "test: " begin
  include("mumps_test.jl")
end
@testset "float complex: " begin
  include("mumps_test_float_complex.jl")
end
@testset " complex: " begin
  include("mumps_test_complex.jl")
end
@testset "mixed: " begin
  include("mumps_test_mixed.jl")
end

MPI.Barrier(comm)
MPI.Finalize()
