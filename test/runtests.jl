using Test
using Random
using LinearAlgebra
using MPI
using MUMPS
using MUMPS: get_sol!
using SparseArrays

Random.seed!(666)  # Random tests are diabolical

@info("MUMPS_INSTALLATION: $(MUMPS.MUMPS_INSTALLATION)")

# include test utilities to obtain quiet ICNTL settings
include("test_utils.jl")

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
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_float.jl")
    end
  end
end
@testset "test: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test.jl")
    end
  end
end
@testset "float complex: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_float_complex.jl")
    end
  end
end
@testset " complex: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_complex.jl")
    end
  end
end
@testset "mixed: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_mixed.jl")
    end
  end
end
@testset "save: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_save.jl")
    end
  end
end
@testset "user permutation: " begin
  redirect_stdout(devnull) do
    redirect_stderr(devnull) do
      include("mumps_test_user_perm.jl")
    end
  end
end

MPI.Barrier(comm)
MPI.Finalize()
