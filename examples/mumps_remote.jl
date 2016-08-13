nprocs() > 1 || error("Run with at least two workers")

@everywhere using MPI
@everywhere using MUMPS
@everywhere MPI.Init()

n = 10
A = rand(n, n)
b = rand(n)
xtrue = A \ b;

# Initialize.
@everywhere mumps64(args...) = Mumps{Float64}(args...)
remote_mumps = remotecall(2, mumps64, mumps_unsymmetric, default_icntl, default_cntl64)

# Associate with matrix A.
MPI.Barrier(MPI.COMM_WORLD)
remotecall(2, associate_matrix!, remote_mumps, A)

# Factorize.
MPI.Barrier(MPI.COMM_WORLD)
remotecall(2, factorize!, remote_mumps)

# Associate with rhs b.
MPI.Barrier(MPI.COMM_WORLD)
remotecall(2, associate_rhs, remote_mumps, b)

# Solve.
MPI.Barrier(MPI.COMM_WORLD)
remotecall(2, solve!, remote_mumps)

# Fetch solution.
MPI.Barrier(MPI.COMM_WORLD)
remote_x = remotecall(2, get_solution, remote_mumps)
x = fetch(remote_x)

MPI.Barrier(MPI.COMM_WORLD)
@printf("Error: %7.1e\n", norm(x - xtrue) / norm(xtrue))

# MPI.Barrier(MPI.COMM_WORLD)
remotecall(2, finalize, fetch(remote_mumps))
