using Base.Test
using MUMPS

# Initialize MPI. Could use MPI.jl here?
mumps_initialize_mpi();

icntl = default_icntl[:];
icntl[1] = 0;
icntl[2] = 0;
icntl[3] = 0;
icntl[4] = 0;

mumps = Mumps(1, icntl);    # Initialize for symmetric definite.

A = spdiagm([1., 2., 3., 4.]);
factorize(mumps, A);  # Analyze and factorize.

rhs = [1., 4., 9., 16.];
x = solve(mumps, rhs);
finalize(mumps);
@test(norm(x - [1, 2, 3, 4]) <= 1.0e-14);

mumps = Mumps(2, icntl);  # General symmetric.

A = rand(4,4); A = sparse(A + A');
factorize(mumps, A);

rhs = rand(4);
x = solve(mumps, rhs);
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);
finalize(mumps);

mumps = Mumps(0, icntl);  # General unsymmetric.

A = sparse(rand(4,4));
factorize(mumps, A);

rhs = rand(4);
x = solve(mumps, rhs);
@test(norm(x - A\rhs)/norm(x) <= 1.0e-12);
finalize(mumps);

mumps_finalize_mpi();

