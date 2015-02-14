icntl = default_icntl[:];
icntl[1] = 0;
icntl[2] = 0;
icntl[3] = 0;
icntl[4] = 0;

mumps1 = Mumps{Float64}(mumps_definite, icntl, default_cntl64);
A = spdiagm([1., 2., 3., 4.]);
factorize(mumps1, A);  # Analyze and factorize.
rhs = [1., 4., 9., 16.];
x = solve(mumps1, rhs);
finalize(mumps1);
MPI.Barrier(comm)
@test(norm(x - [1, 2, 3, 4]) <= 1.0e-14);

mumps2 = Mumps{Float64}(mumps_symmetric, icntl, default_cntl64);
A = rand(4,4); A = sparse(A + A');
factorize(mumps2, A);
rhs = rand(4);
x = solve(mumps2, rhs);
finalize(mumps2);
MPI.Barrier(comm)
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);

mumps3 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64);
A = sparse(rand(4,4));
factorize(mumps3, A);
rhs = rand(4);
x = solve(mumps3, rhs);
finalize(mumps3);
MPI.Barrier(comm)
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);

# Test for solving real div-grad system with single and multiple rhs.
# Based on Lars Ruthotto's initial implementation.

A = get_div_grad(32, 32, 32);
n = size(A, 1);

# Test with single rhs
if MPI.Comm_rank(comm) == root
  println("Test single rhs on div_grad matrix");
end
rhs = randn(n);

x = solve(A, rhs, sym=mumps_definite);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);

# Test with multiple rhs
if MPI.Comm_rank(comm) == root
  println("Test multiple rhs on div_grad matrix");
end
nrhs = 10;
rhs = randn(n, nrhs);

x = solve(A, rhs, sym=mumps_definite);

MPI.Barrier(comm)
relres = zeros(nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]);
  @test(relres[i] <= 1.0e-12);
end

# Test a mixed-type example
A = rand(Float32, 4, 4);
rhs = rand(Int16, 4);
x = solve(A, rhs, sym=mumps_unsymmetric);
MPI.Barrier(comm);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);
