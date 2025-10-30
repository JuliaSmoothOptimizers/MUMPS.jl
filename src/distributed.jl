using MPI

supports_distributed() = Base.find_package("DistributedArrays") !== nothing

function is_darray(A)
  if !supports_distributed()
    return false
  end
  try
    @eval import DistributedArrays
  catch
    return false
  end
  return A isa DistributedArrays.DArray
end

function gather_on_root(A; root::Integer = 0, comm = MPI.COMM_WORLD)
  if !supports_distributed()
    return A
  end
  try
    @eval import DistributedArrays
  catch
    return A
  end

  if A isa DistributedArrays.DArray
    rank = MPI.Comm_rank(comm)
    if rank == root
      return Array(A)
    else
      return nothing
    end
  else
    return A
  end
end
