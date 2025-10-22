# Test mixed-type examples.

# Use deterministic test matrix
A = Float16[1.0 0.5 0.2 0.0; 0.3 2.0 0.5 0.1; 0.0 0.4 3.0 0.5; 0.1 0.0 0.3 4.0]
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
x = solve(A, rhs, sym = mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Here, MUMPS will run in complex single precision.
mumps = quiet_mumps(ComplexF32; sym = mumps_unsymmetric)
x = solve(mumps, A, rhs)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Same, with sparse A.
mumps = quiet_mumps(ComplexF32; sym = mumps_unsymmetric)
x = solve(mumps, sparse(A), rhs)
MPI.Barrier(comm)
finalize(mumps)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))

# Same as above but call associate_matrix! and associate_rhs! directly.
mumps = quiet_mumps(ComplexF32; sym = mumps_unsymmetric)
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
mumps = quiet_mumps(ComplexF32; sym = mumps_unsymmetric)
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
# Use deterministic complex test matrix
A = ComplexF16[1.0+0.1im 0.5 0.2 0.0; 0.3 2.0+0.2im 0.5 0.1; 0.0 0.4 3.0+0.3im 0.5; 0.1 0.0 0.3 4.0+0.4im]
rhs = Array{Int16}([1, 4, 9, 16])
x = solve(A, rhs, sym = mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs)
@test(relres <= eps(Float32)^(1 / 3))
