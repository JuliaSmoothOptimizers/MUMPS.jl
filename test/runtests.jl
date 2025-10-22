using Test
using LinearAlgebra
using MPI
using MUMPS
using MUMPS: get_sol!
using SparseArrays

@info("MUMPS_INSTALLATION: $(MUMPS.MUMPS_INSTALLATION)")

include("get_div_grad.jl")
root = 0

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
@testset "user permutation: " begin
  include("mumps_test_user_perm.jl")
end

MPI.Barrier(comm)
MPI.Finalize()
