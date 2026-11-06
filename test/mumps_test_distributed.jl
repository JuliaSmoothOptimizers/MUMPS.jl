using Test
using LinearAlgebra
using MPI
using MUMPS
using SparseArrays

rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

# keep the entries of A whose position in the list of nonzeros is assigned to this rank
function local_entries(A::SparseMatrixCSC, rank::Integer, nprocs::Integer)
  I, J, V = findnz(A)
  mine = (rank + 1):nprocs:length(V)
  return sparse(I[mine], J[mine], V[mine], size(A)...)
end

function distributed_mumps(::Type{T}, sym) where {T}
  icntl = default_icntl[:]
  icntl[1:4] .= 0
  icntl[18] = 3
  cntl = real(T) == Float32 ? default_cntl32 : default_cntl64
  return Mumps{T}(sym, icntl, cntl)
end

@testset "distributed matrix (ICNTL[18]=3)" begin
  A = get_div_grad(8, 8, 8)
  n = size(A, 1)
  rhs = ones(n)
  tol = sqrt(eps(Float64))

  @testset "$name" for (name, sym, Aref) in (
    ("unsymmetric", mumps_unsymmetric, A + spdiagm(1 => fill(0.5, n - 1))),
    ("symmetric", mumps_symmetric, A),
  )
    mumps = distributed_mumps(Float64, sym)
    A_loc = local_entries(Aref, rank, nprocs)
    associate_matrix!(mumps, A_loc)
    @test mumps.nnz_loc == nnz(sym == mumps_unsymmetric ? A_loc : triu(A_loc))
    factorize!(mumps)
    @test mumps.err == 0
    associate_rhs!(mumps, rhs)
    solve!(mumps)
    @test mumps.err == 0
    x = get_solution(mumps)
    if rank == root
      @test norm(Aref * x - rhs) ≤ tol * norm(rhs) * norm(Aref, 1)
    end
    finalize(mumps)
    MPI.Barrier(comm)
  end

  @testset "duplicate entries are summed" begin
    # every rank holds a share of every entry
    mumps = distributed_mumps(Float64, mumps_symmetric)
    associate_matrix!(mumps, A / nprocs)
    factorize!(mumps)
    associate_rhs!(mumps, rhs)
    solve!(mumps)
    x = get_solution(mumps)
    if rank == root
      @test norm(A * x - rhs) ≤ tol * norm(rhs) * norm(A, 1)
    end
    finalize(mumps)
    MPI.Barrier(comm)
  end

  @testset "unsupported ICNTL[18]" begin
    mumps = distributed_mumps(Float64, mumps_unsymmetric)
    MUMPS.set_icntl!(mumps, 18, 1; displaylevel = 0)
    @test_throws MUMPSException associate_matrix!(mumps, A)
    finalize(mumps)
    MPI.Barrier(comm)
  end
end
