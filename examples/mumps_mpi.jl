using MUMPS
import MPI

MPI.Init()
root = 0;
comm = MPI.COMM_WORLD;

mumps = Mumps(2);  # General symmetric.

# Define problem on the host.
if MPI.Comm_rank(comm) == root
  A = rand(4,4); A = sparse(A + A');
  associate_matrix(mumps, A);
  rhs = rand(4);
  associate_rhs(mumps, rhs);
end

# Solve problem in parallel.
factorize(mumps);
solve(mumps);

MPI.Barrier(comm)

# By default, the solution is assembled and
# overwrites rhs, so only exists on the host.
if MPI.Comm_rank(comm) == root
  x = get_solution(mumps);
  rel_err = norm(x - A\rhs) / norm(x);
  @printf("Error: %7.1e\n", rel_err);
end

finalize(mumps);
MPI.Finalize()
