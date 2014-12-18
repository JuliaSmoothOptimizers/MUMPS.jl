using MUMPS
import MPI

function main(A, rhs)

  MPI.Init()

  mumps = Mumps(2);  # General symmetric.
  factorize(mumps, A);
  x = solve(mumps, rhs);
  finalize(mumps);

  MPI.Barrier(MPI.COMM_WORLD)
  MPI.Finalize()
  return x
end

A = rand(4,4); A = sparse(A + A');
rhs = rand(4);
x = main(A, rhs);
rel_err = norm(x - A\rhs) / norm(x);

# The following usually comes out all garbled.
@printf("Error: %7.1e\n", rel_err);
