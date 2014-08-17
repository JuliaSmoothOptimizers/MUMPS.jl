# Test for solving real div-grad system with single and multiple rhs.
# Based on Lars Ruthotto's initial implementation.

using Base.Test
using MUMPS
include("get_div_grad.jl");

mumps_initialize_mpi();

A = get_div_grad(32, 32, 32);
n = size(A, 1);
icntl = default_icntl[:];
icntl[1] = 0;
icntl[2] = 0;
icntl[3] = 0;
icntl[4] = 0;

# Test with single rhs
println("Test single rhs on div_grad matrix");
rhs = randn(n);

x = solve(A, rhs, 1);
relres = norm(A * x - rhs) / norm(rhs);
@test(relres <= 1.0e-12);

# Test with multiple rhs
println("Test multiple rhs on div_grad matrix");
nrhs = 10;
rhs = randn(n, nrhs);

x = solve(A, rhs, 1);

relres = zeros(nrhs)
for i = 1 : nrhs
  relres[i] =  norm(A * x[:,i] - rhs[:,i]) / norm(rhs[:,i]);
  @test(relres[i] <= 1.0e-12);
end

mumps_finalize_mpi();
