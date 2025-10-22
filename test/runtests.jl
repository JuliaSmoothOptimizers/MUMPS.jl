using Test
using LinearAlgebra
using MPI
using MUMPS
using MUMPS: get_sol!
using SparseArrays

@info("MUMPS_INSTALLATION: $(MUMPS.MUMPS_INSTALLATION)")

# Random tests removed - using deterministic test matrices instead (issue #15)

# Silence helper used in some branches/tests; no-op safe fallback here
if !isdefined(Main, :quiet_mumps)
  function quiet_mumps(f::Function)
    Base.redirect_stdout(devnull) do
      Base.redirect_stderr(devnull) do
        f()
      end
    end
  end
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
@testset "save: " begin
  include("mumps_test_save.jl")
end
@testset "user permutation: " begin
  include("mumps_test_user_perm.jl")
end

MPI.Barrier(comm)
MPI.Finalize()
