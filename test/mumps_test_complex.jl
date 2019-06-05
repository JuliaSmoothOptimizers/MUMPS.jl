icntl = default_icntl[:]
icntl[1] = 0
icntl[2] = 0
icntl[3] = 0
icntl[4] = 0
tol = sqrt(eps(Float64))

mumps1 = Mumps{ComplexF64}(mumps_definite, icntl, default_cntl64)
A = sparse(Diagonal(([1., 2., 3., 4.])))
factorize!(mumps1, A)  # Analyze and factorize.
rhs = [1., 4., 9., 16.]
x = solve(mumps1, rhs)
finalize(mumps1)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

mumps2 = Mumps{ComplexF64}(mumps_symmetric, icntl, default_cntl64)
A = random_matrix(Float64, [1, 2, 3, 4], 4, 4); A = sparse(A + A')
factorize!(mumps2, A)
rhs = [1., 4., 9., 16.]
x = solve(mumps2, rhs)
finalize(mumps2)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

mumps3 = Mumps{ComplexF64}(mumps_unsymmetric, icntl, default_cntl64)
A = sparse(random_matrix(ComplexF64, [1, 2, 3, 4], 4, 4))
factorize!(mumps3, A)
rhs = [1., 4., 9., 16.] + im * [1., 4., 9., 16.]
x = solve(mumps3, rhs)
finalize(mumps3)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

# Test convenience interface.

n = 10
n3 = n * n * n
A = convert(SparseMatrixCSC{ComplexF64,Int32}, map(ComplexF64, get_div_grad(n, n, n)))

# Test with single rhs
if MPI.Comm_rank(comm) == root
  println("Test single rhs on div_grad matrix")
end
rhs = ones(n3) + im * ones(n3)

x = solve(A, rhs, sym=mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs) / norm(A, 1)
@test(relres <= tol);

# Test with multiple rhs
if MPI.Comm_rank(comm) == root
  println("Test multiple rhs on div_grad matrix")
end
nrhs = 5;
rhs = map(ComplexF64, ones(n3, nrhs) + im * ones(n3, nrhs)) * diagm(0 => Array{Float64}(collect(1:nrhs)))

x = solve(A, rhs, sym=mumps_unsymmetric)

MPI.Barrier(comm)
relres = zeros(Float64, nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]) /  norm(A, 1)
  @test(relres[i] <= 1.0e-8)
end
