# Test mixed-type examples.

# For now both A and rhs will be converted to Float64.
A = random_matrix(Float16, [1, 2, 3, 4], 4, 4);
rhs = Array{Float32}([1., 4., 9., 16.]);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);

# Here, MUMPS will run in complex single precision.
mumps = Mumps{Complex64}(mumps_unsymmetric, default_icntl, default_cntl32);
x = solve(mumps, A, rhs);
MPI.Barrier(comm);
finalize(mumps);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-6);

# Same, with sparse A.
mumps = Mumps{Complex64}(mumps_unsymmetric, default_icntl, default_cntl32);
x = solve(mumps, sparse(A), rhs);
MPI.Barrier(comm);
finalize(mumps);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-6);

# Same as above but call associate_matrix! and associate_rhs! directly.
mumps = Mumps{Complex64}(mumps_unsymmetric, default_icntl, default_cntl32);
associate_matrix!(mumps, A);
factorize!(mumps);
associate_rhs!(mumps, rhs);
solve!(mumps);
x = get_solution(mumps);
MPI.Barrier(comm);
finalize(mumps);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-6);

# And the same with sparse A.
mumps = Mumps{Complex64}(mumps_unsymmetric, default_icntl, default_cntl32);
associate_matrix!(mumps, sparse(A));
factorize!(mumps);
associate_rhs!(mumps, rhs);
solve!(mumps);
x = get_solution(mumps);
MPI.Barrier(comm);
finalize(mumps);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-6);

# For now both A and rhs will be converted to Complex128.
A = random_matrix(Complex32, [1, 2, 3, 4], 4, 4);
rhs = Array{Int16}([1, 4, 9, 16]);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);
