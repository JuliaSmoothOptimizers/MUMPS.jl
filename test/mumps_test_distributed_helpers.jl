using Test
using MPI
using MUMPS

if Base.find_package("DistributedArrays") === nothing
  @info "DistributedArrays not installed; skipping distributed helper tests"
  @testset "distributed helpers: skipped (no DistributedArrays)" begin
    @test true
  end
else
  import DistributedArrays

  comm = MPI.COMM_WORLD
  rank = MPI.Comm_rank(comm)

  @testset "distributed helpers" begin
    localA = reshape(collect(1:12), 3, 4)
    darr = DistributedArrays.distribute(localA)

    g = MUMPS.gather_on_root(darr; root = 0, comm = comm)

    if rank == 0
      @test isa(g, Array)
      @test size(g) == size(localA)
      @test all(g .== localA)
    else
      @test g === nothing
    end
  end
end
