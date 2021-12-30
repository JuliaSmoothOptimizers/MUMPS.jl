using Test
using LinearAlgebra
using MPI
using MUMPS
using SparseArrays

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

# 1D finite difference on staggered grid
ddx(n :: Int) = sparse(Bidiagonal(-ones(n), ones(n-1), :U))

# Based on Lars Ruthotto's initial implementation.
function get_div_grad(n1 :: Int, n2 :: Int, n3 :: Int)

  # Divergence
  In1 = sparse(1.0I, n1, n1)
  In2 = sparse(1.0I, n2, n2)
  In3 = sparse(1.0I, n3, n3)
  D1 = kron(In3, kron(In2, ddx(n1)))
  D2 = kron(In3, kron(ddx(n2),   In1))
  D3 = kron(ddx(n3),   kron(In2, In1))

  # DIV from faces to cell-centers
  Div = [D1 D2 D3]

  return Div * Div'
end

# Initialize MPI.
MPI.Init()
comm = MPI.COMM_WORLD
root = 0

@testset "float: " begin include("mumps_test_float.jl") end
@testset "test: " begin include("mumps_test.jl") end
@testset "float complex: " begin include("mumps_test_float_complex.jl") end
@testset " complex: " begin include("mumps_test_complex.jl") end
@testset "mixed: " begin include("mumps_test_mixed.jl") end

MPI.Barrier(comm)
MPI.Finalize()
