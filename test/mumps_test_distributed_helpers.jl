using Test
using MPI
using MUMPS
import DistributedArrays

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

@testset "distributed helpers" begin
  localA = reshape(collect(1:12), 3, 4)
  darr = DistributedArrays.distribute(localA)

  if rank == 0
    g = Array(darr)
  else
    g = nothing
  end

  if rank == 0
    @test isa(g, Array)
    @test size(g) == size(localA)
    @test all(g .== localA)
  else
    @test g === nothing
  end
end
