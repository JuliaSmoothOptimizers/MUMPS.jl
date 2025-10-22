using Test
using LinearAlgebra
using MPI
using MUMPS
using SparseArrays

@testset "ScaLAPACK enabled (ICNTL[13]=0)" begin
  n = 20
  A = sprand(Float64, n, n, 0.3) + I
  A = A + A'
  rhs = rand(Float64, n)

  icntl = default_icntl[:]
  icntl[1:4] .= 0
  icntl[13] = 0

  mumps = Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)

  @test mumps.icntl[13] == 0

  factorize!(mumps, A)
  x = solve(mumps, rhs)

  tol = sqrt(eps(Float64))
  @test norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1)

  finalize(mumps)
  MPI.Barrier(comm)
end

@testset "ScaLAPACK disabled (ICNTL[13]>0)" begin
  n = 20
  A = sprand(Float64, n, n, 0.3) + I
  A = A + A'
  rhs = rand(Float64, n)

  icntl = default_icntl[:]
  icntl[1:4] .= 0
  icntl[13] = 999

  mumps = Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)

  @test mumps.icntl[13] == 999

  factorize!(mumps, A)
  x = solve(mumps, rhs)

  tol = sqrt(eps(Float64))
  @test norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1)

  finalize(mumps)
  MPI.Barrier(comm)
end

@testset "ScaLAPACK with unsymmetric matrix" begin
  n = 20
  A = sprand(Float64, n, n, 0.3) + I
  rhs = rand(Float64, n)

  icntl = default_icntl[:]
  icntl[1:4] .= 0
  icntl[13] = 0

  mumps = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)

  @test mumps.icntl[13] == 0

  factorize!(mumps, A)
  x = solve(mumps, rhs)

  tol = sqrt(eps(Float64))
  @test norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1)

  finalize(mumps)
  MPI.Barrier(comm)
end

@testset "Modify ICNTL[13] with set_icntl!" begin
  n = 20
  A = sprand(Float64, n, n, 0.3) + I
  A = A + A'
  rhs = rand(Float64, n)

  icntl = default_icntl[:]
  icntl[1:4] .= 0

  mumps = Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)

  MUMPS.set_icntl!(mumps, 13, 10; displaylevel = 0)
  @test mumps.icntl[13] == 10

  factorize!(mumps, A)
  x = solve(mumps, rhs)

  tol = sqrt(eps(Float64))
  @test norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1)

  finalize(mumps)
  MPI.Barrier(comm)
end
