# Test mixed-type examples.

A = random_matrix(Float16, [1, 2, 3, 4], 4, 4)
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
x = solve(A, rhs, sym = mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Here, MUMPS will run in complex single precision.
mumps = Mumps{ComplexF32}(mumps_unsymmetric, default_icntl, default_cntl32)
x = solve(mumps, A, rhs)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Same, with sparse A.
mumps = Mumps{ComplexF32}(mumps_unsymmetric, default_icntl, default_cntl32)
x = solve(mumps, sparse(A), rhs)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Same as above but call associate_matrix! and associate_rhs! directly.
mumps = Mumps{ComplexF32}(mumps_unsymmetric, default_icntl, default_cntl32)
associate_matrix!(mumps, A)
factorize!(mumps)
associate_rhs!(mumps, rhs)
solve!(mumps)
x = get_solution(mumps)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# And the same with sparse A.
mumps = Mumps{ComplexF32}(mumps_unsymmetric, default_icntl, default_cntl32)
associate_matrix!(mumps, sparse(A))
factorize!(mumps)
associate_rhs!(mumps, rhs)
solve!(mumps)
x = get_solution(mumps)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# For now both A and rhs will be converted to ComplexF64.
A = random_matrix(ComplexF16, [1, 2, 3, 4], 4, 4)
rhs = Array{Int16}([1, 4, 9, 16])
x = solve(A, rhs, sym = mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))
