"""
Enumeration of MUMPS job types; see MUMPS documentation for details.
"""
@enum MUMPS_JOB::MUMPS_INT begin
  INITIALIZE = -1 # initializes an instance of MUMPS
  TERMINATE = -2 # terminates an instance of MUMPS
  DELETE_DATA = -3 # removes data saved to disk
  FACTOR_CLEANUP = -4 # frees all internal data, except those from analysis
  # SURPRESS = -200 # "(experimental, subject to change) surpresses all MUMPS 
  # out-of-core factor files associated with MPI processes and returns"
  ANALYZE = 1 # performs the analysis phase
  FACTOR = 2 # performs the factorization phase
  SOLVE = 3 # performs the solve phase
  ANALYZE_FACTOR = 4 # combines analyze and factorization phases
  FACTOR_SOLVE = 5 # combines factorization and solve phases
  ANALYZE_FACTOR_SOLVE = 6 # combines analyze, factorization, and solve phases
  SAVE_DATA = 7 # save/restore feature; save internal mumps data to disk.
  RESTORE_DATA = 8 # save/restore feature; restore internal mumps data from disk.
  # RHS_DISTRIBUTION = 9 # "computes before the solution phase a possible distribution
  # for the right-hand sides [among the processors]"
end

const SOLVE_JOBS = (SOLVE, FACTOR_SOLVE, ANALYZE_FACTOR_SOLVE)
const ONLY_FACTORED = (FACTOR, ANALYZE_FACTOR)

is_factored(job::MUMPS_JOB) = job >= FACTOR

function get_name(job::MUMPS_JOB)
  if job == INITIALIZE
    return "initialize"
  elseif job == TERMINATE
    return "terminate"
  elseif job == DELETE_DATA
    return "delete saved data"
  elseif job == FACTOR_CLEANUP
    return "factorization cleanup [free all internal data except analysis]"
  elseif job == ANALYZE
    return "analyze"
  elseif job == FACTOR
    return "factorize"
  elseif job == SOLVE
    return "solve"
  elseif job == ANALYZE_FACTOR
    return "analyze + factorize"
  elseif job == FACTOR_SOLVE
    return "factorize + solve"
  elseif job == ANALYZE_FACTOR_SOLVE
    return "analyze + factorize + solve"
  elseif job == SAVE_DATA
    return "save"
  elseif job == RESTORE_DATA
    return "restore"
  else
    return "unrecognized"
  end
end
