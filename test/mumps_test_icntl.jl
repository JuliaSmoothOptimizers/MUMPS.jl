# Tests for the control parameters introduced up to MUMPS 5.9.

tol = sqrt(eps(Float64))
A = sparse([4.0 1.0 0.0 0.0; 1.0 4.0 1.0 0.0; 0.0 1.0 4.0 1.0; 0.0 0.0 1.0 4.0])
rhs = [1.0, 2.0, 3.0, 4.0]

@testset "default control parameters" begin
  @test length(default_icntl) == 60
  @test default_icntl[36] == 1
  @test default_icntl[38] == 600
  @test default_icntl[48] == 1
  @test default_icntl[58] == 2

  mumps = quiet_mumps(Float64)
  MUMPS.set_icntl!(mumps, 48, 0; displaylevel = 0)
  MUMPS.default_icntl!(mumps)
  @test collect(mumps.icntl) == default_icntl
  finalize(mumps)

  # CNTL(1) = -1 (automatic) is forwarded to MUMPS and the defaults are not mutated.
  mumps = Mumps{Float64}(mumps_definite, quiet_icntl(), default_cntl64)
  @test mumps.cntl[1] == -1
  @test default_cntl64[1] == -1
  finalize(mumps)
  MPI.Barrier(comm)
end

@testset "get_icntl keywords" begin
  icntl = get_icntl(
    blr = 1,
    nthreads = 2,
    tree_parallelism = false,
    rank_revealing = true,
    single_precision_factorization = true,
    null_pivots = true,
  )
  @test length(icntl) == 60
  @test icntl[35] == 1
  @test icntl[16] == 2
  @test icntl[48] == 0
  @test icntl[56] == 1
  @test icntl[47] == 1
  @test icntl[24] == 1
  @test get_icntl()[48] == 1
  @test_throws ArgumentError get_icntl(blr = 4)
  @test_throws ArgumentError get_icntl(nthreads = -1)
end

@testset "display" begin
  io = IOBuffer()
  MUMPS.display_icntl(io, default_icntl)
  out = String(take!(io))
  @test length(split(strip(out), '\n')) == 60
  @test occursin("rank-revealing", out)
  @test occursin("tree parallelism", out)
  @test occursin("single precision factorization", out)

  MUMPS.display_cntl(io, default_cntl64)
  out = String(take!(io))
  @test occursin("automatic", out)

  mumps = quiet_mumps(Float64)
  MUMPS.display_cntl(io, mumps)
  @test !isempty(String(take!(io)))
  show(io, mumps)
  @test !isempty(String(take!(io)))
  finalize(mumps)
  MPI.Barrier(comm)
end

@testset "solve without tree parallelism" begin
  mumps = Mumps{Float64}(mumps_symmetric, get_icntl(tree_parallelism = false), default_cntl64)
  @test mumps.icntl[48] == 0
  x = solve(mumps, A, rhs)
  finalize(mumps)
  MPI.Barrier(comm)
  @test norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1)
end

@testset "single precision factorization (ICNTL(47))" begin
  mumps = Mumps{Float64}(mumps_unsymmetric, get_icntl(itref = 3), default_cntl64)
  MUMPS.single_precision_factorization!(mumps)
  @test mumps.icntl[47] == 1
  x = solve(mumps, A, rhs)
  finalize(mumps)
  MPI.Barrier(comm)
  @test norm(A * x - rhs) <= sqrt(eps(Float32)) * norm(rhs) * norm(A, 1)

  mumps32 = quiet_mumps(Float32)
  @test_throws MUMPSException MUMPS.single_precision_factorization!(mumps32)
  finalize(mumps32)
  MPI.Barrier(comm)
end
