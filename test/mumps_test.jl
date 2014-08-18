using Base.Test
using MUMPS

include("get_div_grad.jl");

# Initialize MPI. Could also use MPI.jl.
mumps_initialize_mpi();

icntl = default_icntl[:];
#=icntl[1] = 0;=#
#=icntl[2] = 0;=#
#=icntl[3] = 0;=#
#=icntl[4] = 0;=#

mumps1 = Mumps(1, icntl=icntl);    # Initialize for symmetric definite.
A = spdiagm([1., 2., 3., 4.]);
factorize(mumps1, A);  # Analyze and factorize.
rhs = [1., 4., 9., 16.];
x = solve(mumps1, rhs);
finalize(mumps1);
@test(norm(x - [1, 2, 3, 4]) <= 1.0e-14);

mumps2 = Mumps(2, icntl=icntl);  # General symmetric.
A = rand(4,4); A = sparse(A + A');
factorize(mumps2, A);
rhs = rand(4);
x = solve(mumps2, rhs);
finalize(mumps2);
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);

mumps3 = Mumps(0, icntl=icntl);  # General unsymmetric.
A = sparse(rand(4,4));
factorize(mumps3, A);
rhs = rand(4);
x = solve(mumps3, rhs);
finalize(mumps3);
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);

# Test for solving real div-grad system with single and multiple rhs.
# Based on Lars Ruthotto's initial implementation.

A = get_div_grad(32, 32, 32);
n = size(A, 1);

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

