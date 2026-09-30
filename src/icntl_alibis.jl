# this file contains a bunch of shortcuts for manipulating ICNTL.

# export set_error_stream!, set_diagnostics_stream!, set_info_stream!, set_print_level!,
# suppress_printing!, toggle_printing!, suppress_display!, toggle_display!,
# sparse_matrix!, dense_matrix!,
# sparse_rhs!, dense_rhs!,
# toggle_null_pivot!

const ICNTL_DEFAULT = (
  6,
  0,
  6,
  2,
  0,
  7,
  7,
  77,
  1,
  0,
  0,
  0,
  0,
  20,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  -32,
  1,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  1,
  0,
  600,
  500,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  1,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  2,
  0,
  0,
)
"""
    default_icntl!(mumps)

reset ICNTL to its default
"""
function default_icntl!(mumps::Mumps)
  mumps.icntl = ICNTL_DEFAULT
  return nothing
end

set_error_stream!(mumps::Mumps, i) = begin
  set_icntl!(mumps, 1, i; displaylevel = 0)
  return mumps
end
set_diagnostics_stream!(mumps::Mumps, i) = begin
  set_icntl!(mumps, 2, i; displaylevel = 0)
  return mumps
end
set_info_stream!(mumps::Mumps, i) = begin
  set_icntl!(mumps, 3, i; displaylevel = 0)
  return mumps
end
set_print_level!(mumps::Mumps, i) = begin
  set_icntl!(mumps, 4, i; displaylevel = 0)
  return mumps
end

suppress_printing!(mumps::Mumps) = begin
  set_print_level!(mumps, 1)
  return mumps
end
toggle_printing!(mumps::Mumps) = begin
  set_print_level!(mumps, mod1(mumps.icntl[4] + 1, 2))
  return mumps
end
suppress_display! = suppress_printing!
toggle_display! = toggle_printing!

sparse_matrix!(mumps::Mumps) = begin
  set_icntl!(mumps, 5, 0; displaylevel = 0)
  return mumps
end
dense_matrix!(mumps::Mumps) = begin
  set_icntl!(mumps, 5, 1; displaylevel = 0)
  return mumps
end

LinearAlgebra.transpose!(mumps::Mumps) = begin
  set_icntl!(mumps, 9, mod(mumps.icntl[9] + 1, 2); displaylevel = 0)
  return mumps
end

sparse_rhs!(mumps::Mumps) = begin
  set_icntl!(mumps, 20, 1; displaylevel = 0)
  return mumps
end
dense_rhs!(mumps::Mumps) = begin
  set_icntl!(mumps, 20, 0; displaylevel = 0)
  return mumps
end

toggle_null_pivot!(mumps::Mumps) = begin
  set_icntl!(mumps, 24, mod(mumps.icntl[24] + 1, 2); displaylevel = 0)
  return mumps
end

set_num_threads!(mumps::Mumps, n::Integer) = begin
  set_icntl!(mumps, 16, n; displaylevel = 0)
  return mumps
end

toggle_tree_parallelism!(mumps::Mumps) = begin
  set_icntl!(mumps, 48, mod(mumps.icntl[48] + 1, 2); displaylevel = 0)
  return mumps
end

toggle_rank_revealing!(mumps::Mumps) = begin
  set_icntl!(mumps, 56, mod(mumps.icntl[56] + 1, 2); displaylevel = 0)
  return mumps
end

"""
    set_blr!(mumps, tol; mode=1)

Activate the Block Low-Rank factorization (ICNTL(35)=`mode`) with dropping parameter
CNTL(7)=`tol`. Use `mode=0` to deactivate it.
"""
function set_blr!(mumps::Mumps, tol::AbstractFloat; mode::Integer = 1)
  0 ≤ mode ≤ 3 || throw(MUMPSException("ICNTL(35) must be 0, 1, 2 or 3, got $mode"))
  set_icntl!(mumps, 35, mode; displaylevel = 0)
  set_cntl!(mumps, 7, tol; displaylevel = 0)
  return mumps
end

"""
    single_precision_factorization!(mumps, flag=true)

Perform the factorization in single precision within a double precision instance
(ICNTL(47), MUMPS ≥ 5.9). Iterative refinement (ICNTL(10)) is recommended to recover
double precision accuracy.
"""
function single_precision_factorization!(mumps::Mumps{TC, TR}, flag::Bool = true) where {TC, TR}
  TR == Float64 ||
    throw(MUMPSException("ICNTL(47) requires a Float64 or ComplexF64 instance, got $TC"))
  set_icntl!(mumps, 47, flag ? 1 : 0; displaylevel = 0)
  return mumps
end
