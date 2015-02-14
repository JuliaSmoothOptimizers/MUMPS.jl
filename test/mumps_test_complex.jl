icntl = default_icntl[:];
icntl[1] = 0;
icntl[2] = 0;
icntl[3] = 0;
icntl[4] = 0;

mumps1 = Mumps{Complex128}(mumps_definite, icntl, default_cntl64);
A = spdiagm([1., 2., 3., 4.]);
factorize(mumps1, A);  # Analyze and factorize.
rhs = [1., 4., 9., 16.];
x = solve(mumps1, rhs);
finalize(mumps1);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= 1.0e-12 * norm(rhs));

mumps2 = Mumps{Complex128}(mumps_symmetric, icntl, default_cntl64);
A = rand(4,4); A = sparse(A + A');
factorize(mumps2, A);
rhs = rand(4);
x = solve(mumps2, rhs);
finalize(mumps2);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= 1.0e-12 * norm(rhs));

mumps3 = Mumps{Complex128}(mumps_unsymmetric, icntl, default_cntl64);
A = sparse(rand(4,4)) + im * sparse(rand(4,4));
factorize(mumps3, A);
rhs = rand(4) + im * rand(4);
x = solve(mumps3, rhs);
finalize(mumps3);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= 1.0e-12 * norm(rhs));

# Test convenience interface.

n = 100;
A = sprand(100, 100, .2) + im * sprand(100, 100, .2);

# Test with single rhs
if MPI.Comm_rank(comm) == root
  println("Test single rhs on div_grad matrix");
end
rhs = randn(n) + im * randn(n);

x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);

# Test with multiple rhs
if MPI.Comm_rank(comm) == root
  println("Test multiple rhs on div_grad matrix");
end
nrhs = 10;
rhs = randn(n, nrhs) + im * randn(n, nrhs);

x = solve(A, rhs, sym=mumps_unsymmetric);

MPI.Barrier(comm)
relres = zeros(Float64, nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]);
  @test(relres[i] <= 1.0e-8);
end

# Test a mixed-type example
A = rand(Complex32, 4, 4);
rhs = rand(Int16, 4);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);
