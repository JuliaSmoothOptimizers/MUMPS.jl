icntl = default_icntl[:];
icntl[1] = 0;
icntl[2] = 0;
icntl[3] = 0;
icntl[4] = 0;
tol = sqrt(eps(Float64))

mumps1 = Mumps{Float64}(mumps_definite, icntl, default_cntl64);
A = sparse(Diagonal([1., 2., 3., 4.]))
factorize!(mumps1, A);  # Analyze and factorize.
rhs = [1., 4., 9., 16.];
x = solve(mumps1, rhs);
finalize(mumps1);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1));

mumps2 = Mumps{Float64}(mumps_symmetric, icntl, default_cntl64);
A = random_matrix(Float64, [1, 2, 3, 4], 4, 4); A = sparse(A + A');
factorize!(mumps2, A);
rhs = [1., 4., 9., 16.];
x = solve(mumps2, rhs);
finalize(mumps2);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1));

mumps3 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64);
A = sparse(random_matrix(Float64, [1, 2, 3, 4], 4, 4));
factorize!(mumps3, A);
rhs = [1., 4., 9., 16.];
x = solve(mumps3, rhs);
finalize(mumps3);
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1));

# Test for solving real div-grad system with single and multiple rhs.
# Based on Lars Ruthotto's initial implementation.

n = 10
n3 = 10 * 10 * 10
A = convert(SparseMatrixCSC{Float64,Int32}, get_div_grad(n, n, n))

# Test with single rhs
if MPI.Comm_rank(comm) == root
  println("Test single rhs on div_grad matrix");
end
rhs = ones(n3);

x = solve(A, rhs, sym=mumps_definite);
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= tol);

# Test with multiple rhs
if MPI.Comm_rank(comm) == root
  println("Test multiple rhs on div_grad matrix");
end
nrhs = 5;
rhs = ones(n3, nrhs) * diagm(0 => collect(1:nrhs))

x = solve(A, rhs, sym=mumps_definite);

MPI.Barrier(comm)
relres = zeros(nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]) / norm(A, 1);
  @test(relres[i] <= tol);
end
