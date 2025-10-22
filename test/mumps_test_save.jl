# Use deterministic sparse symmetric test matrix instead of random
A = sparse([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 3, 4, 5, 6, 7, 8, 9, 10],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9],
           [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
           10, 10)
S = A + A'
mumps = quiet_mumps(Float64; sym = mumps_symmetric)
associate_matrix!(mumps, S)
MUMPS.set_job!(mumps, 1)
MUMPS.invoke_mumps!(mumps)
save_dir = mktempdir(; prefix = "_mumps_test_save_")
MUMPS.set_save_dir!(mumps, save_dir)
MUMPS.set_job!(mumps, 7)
MUMPS.invoke_mumps!(mumps)
@test length(readdir(save_dir)) > 0
MUMPS.set_job!(mumps, -3)
MUMPS.invoke_mumps!(mumps)
finalize(mumps)
