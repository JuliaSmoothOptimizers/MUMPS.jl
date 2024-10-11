icntl = get_icntl(det = true, ooc = true, itref = 1);
tol = sqrt(eps(Float32))

mumps1 = Mumps{Float32}(mumps_definite, icntl, default_cntl32);
A = sparse(Diagonal(Array{Float32}([1.0, 2.0, 3.0, 4.0])))
factorize!(mumps1, A);  # Analyze and factorize.
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
x = solve(mumps1, rhs)
finalize(mumps1)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

mumps1_unsafe = Mumps{Float32}(mumps_definite, icntl, default_cntl32);
A = sparse(Diagonal(Array{Float32}([1.0, 2.0, 3.0, 4.0])))
associate_matrix!(mumps1_unsafe, A; unsafe=true)
factorize!(mumps1_unsafe);  # Analyze and factorize.
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
orig_rhs = copy(rhs)
associate_rhs!(mumps1_unsafe, rhs; unsafe=true)
solve!(mumps1_unsafe)
x = similar(orig_rhs)
get_sol!(x, mumps1_unsafe)
finalize(mumps1_unsafe)
MPI.Barrier(comm)
@test(norm(A * x - orig_rhs) <= tol * norm(orig_rhs) * norm(A, 1))

mumps2 = Mumps{Float32}(mumps_symmetric, icntl, default_cntl32)
A = random_matrix(Float32, [1, 2, 3, 4], 4, 4);
A = sparse(A + A');
factorize!(mumps2, A)
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
x = solve(mumps2, rhs)
finalize(mumps2)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

mumps3 = Mumps{Float32}(mumps_unsymmetric, icntl, default_cntl32)
A = sparse(random_matrix(Float32, [1, 2, 3, 4], 4, 4))
factorize!(mumps3, A)
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
x = solve(mumps3, rhs)
finalize(mumps3)
MPI.Barrier(comm)
@test(norm(A * x - rhs) <= tol * norm(rhs) * norm(A, 1))

mumps3_unsafe = Mumps{Float32}(mumps_unsymmetric, icntl, default_cntl32);
A = sparse(random_matrix(Float32, [1, 2, 3, 4], 4, 4))
associate_matrix!(mumps3_unsafe, A; unsafe=true)
factorize!(mumps3_unsafe);  # Analyze and factorize.
rhs = Array{Float32}([1.0, 4.0, 9.0, 16.0])
orig_rhs = copy(rhs)
associate_rhs!(mumps3_unsafe, rhs; unsafe=true)
solve!(mumps3_unsafe)
x = similar(orig_rhs)
get_sol!(x, mumps3_unsafe)
finalize(mumps3_unsafe)
MPI.Barrier(comm)
@test(norm(A * x - orig_rhs) <= tol * norm(orig_rhs) * norm(A, 1))

# Test convenience interface.

n = 10
n3 = n * n * n
A = convert(SparseMatrixCSC{Float64, Int32}, get_div_grad(n, n, n))

# Test with single rhs
if MPI.Comm_rank(comm) == root
  println("Test single rhs")
end
rhs = ones(n3);

# For now everything will be converted to Float64.
x = solve(A, rhs, sym = mumps_unsymmetric)
MPI.Barrier(comm)
relres = norm(A * x - rhs) / norm(rhs) / norm(A, 1)
@test(relres <= 1.0e-12)

# Test with multiple rhs
if MPI.Comm_rank(comm) == root
  println("Test multiple rhs")
end
nrhs = 5;
rhs = ones(n3, nrhs) * diagm(0 => collect(1:nrhs))

x = solve(A, rhs, sym = mumps_unsymmetric);

MPI.Barrier(comm)
relres = zeros(Float32, nrhs)
for i = 1:nrhs
  relres[i] = norm(A * x[:, i] - rhs[:, i]) / norm(rhs[:, i]) / norm(A, 1)
  @test(relres[i] <= tol)
end
