# Test for solving real div-grad system with single and multiple rhs.
# Based on Lars Ruthotto's initial implementation.

using MUMPS
include("get_div_grad.jl");

mumps_initialize_mpi();

A = get_div_grad(32, 32, 32);
n = size(A, 1);

# Test with single rhs
println("Test for real SPD matrix: one rhs");
rhs = randn(n);

x = solve(A, rhs, 1);
relres = norm(A * x - rhs) / norm(rhs);
@printf("Relative residual: %7.1e\n", relres);

# Test with multiple rhs
println("Test for real SPD matrix: multiple rhs");
nrhs = 10;
rhs = randn(n, nrhs);

x = solve(A, rhs, 1);

relres = zeros(nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]);
  @printf("Relative residual %2d: %7.1e\n", i, relres[i]);
end

mumps_finalize_mpi();
