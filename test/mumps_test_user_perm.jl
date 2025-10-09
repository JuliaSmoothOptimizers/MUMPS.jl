# Test user-supplied permutation functionality

icntl = default_icntl[:]
icntl[1] = 0  # No error output
icntl[2] = 0  # No diagnostic output  
icntl[3] = 0  # No global info output
icntl[4] = 0  # No output level
tol = sqrt(eps(Float64))

# Test with a simple diagonal matrix where permutation effects are clear
@testset "User-supplied permutation" begin
  # Create a simple test matrix
  A = sparse([
    1.0 0.0 0.0 0.0;
    0.0 2.0 0.0 0.0;
    0.0 0.0 3.0 0.0;
    0.0 0.0 0.0 4.0
  ])
  rhs = [1.0, 4.0, 9.0, 16.0]
  expected_x = [1.0, 2.0, 3.0, 4.0]

  # Test 1: Default ordering (should work)
  mumps1 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  x1 = solve(mumps1, A, rhs)
  finalize(mumps1)
  MPI.Barrier(comm)
  @test(norm(A * x1 - rhs) <= tol * norm(rhs) * norm(A, 1))

  # Test 2: Identity permutation (should give same result)
  mumps2 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  identity_perm = [1, 2, 3, 4]
  set_user_perm!(mumps2, identity_perm)
  # Verify ICNTL[7] was set correctly
  @test mumps2.icntl[7] == 1
  x2 = solve(mumps2, A, rhs)
  finalize(mumps2)
  MPI.Barrier(comm)
  @test(norm(A * x2 - rhs) <= tol * norm(rhs) * norm(A, 1))
  @test(norm(x2 - expected_x) <= tol * norm(expected_x))

  # Test 3: Reverse permutation
  mumps3 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  reverse_perm = [4, 3, 2, 1]
  set_user_perm!(mumps3, reverse_perm)
  @test mumps3.icntl[7] == 1
  x3 = solve(mumps3, A, rhs)
  finalize(mumps3)
  MPI.Barrier(comm)
  @test(norm(A * x3 - rhs) <= tol * norm(rhs) * norm(A, 1))

  # Test 4: Test with get_icntl convenience function
  icntl_user = get_icntl(user_perm = true)
  @test icntl_user[7] == 1
  mumps4 = Mumps{Float64}(mumps_unsymmetric, icntl_user, default_cntl64)
  set_user_perm!(mumps4, identity_perm)
  x4 = solve(mumps4, A, rhs)
  finalize(mumps4)
  MPI.Barrier(comm)
  @test(norm(A * x4 - rhs) <= tol * norm(rhs) * norm(A, 1))

  # Test 5: Test unsafe option
  mumps5 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  perm_unsafe = [1, 2, 3, 4]
  set_user_perm!(mumps5, perm_unsafe; unsafe = true)
  @test mumps5.icntl[7] == 1
  x5 = solve(mumps5, A, rhs)
  finalize(mumps5)
  MPI.Barrier(comm)
  @test(norm(A * x5 - rhs) <= tol * norm(rhs) * norm(A, 1))

  # Test 6: Test with complex matrix
  A_complex = sparse(complex.([
    1.0 0.0 0.0 0.0;
    0.0 2.0 0.0 0.0;
    0.0 0.0 3.0 0.0;
    0.0 0.0 0.0 4.0
  ]))
  rhs_complex = complex([1.0, 4.0, 9.0, 16.0])

  mumps6 = Mumps{ComplexF64}(mumps_unsymmetric, icntl, default_cntl64)
  set_user_perm!(mumps6, identity_perm)
  @test mumps6.icntl[7] == 1
  x6 = solve(mumps6, A_complex, rhs_complex)
  finalize(mumps6)
  MPI.Barrier(comm)
  @test(norm(A_complex * x6 - rhs_complex) <= tol * norm(rhs_complex) * norm(A_complex, 1))

  # Test 7: Test error handling for wrong permutation size
  mumps7 = Mumps{Float64}(mumps_unsymmetric, icntl, default_cntl64)
  wrong_size_perm = [1, 2, 3]  # Too short for 4x4 matrix
  set_user_perm!(mumps7, wrong_size_perm)  # This should still work, MUMPS will handle the error
  # We don't test the solve here since it would fail, just test that the function completes
  @test mumps7.icntl[7] == 1
  finalize(mumps7)
  MPI.Barrier(comm)
end
