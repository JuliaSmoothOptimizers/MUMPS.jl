# Test mixed-type examples.

# For now both A and rhs will be converted to Float64.
A = rand(Float16, 4, 4);
rhs = rand(Int32, 4);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);

# For now both A and rhs will be converted to Complex128.
A = rand(Complex32, 4, 4);
rhs = rand(Int16, 4);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);
