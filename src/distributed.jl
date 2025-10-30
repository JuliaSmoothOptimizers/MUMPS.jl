"""Helpers for optional support of DistributedArrays/DArray without adding a
hard dependency. These helpers detect `DistributedArrays` at runtime and provide
light-weight utilities such as `gather_on_root` which will gather a distributed
array onto the specified root rank (using `MPI.COMM_WORLD` by default).

This file intentionally avoids adding `DistributedArrays` to `Project.toml`.
If `DistributedArrays` is available at runtime the functions below will use it,
otherwise they fall back to no-op behaviour.
"""

using MPI

"""Return true if `DistributedArrays` is available in the current environment."""
supports_distributed() = Base.find_package("DistributedArrays") !== nothing

"""Return true if `A` is a DistributedArrays.DArray.

If `DistributedArrays` is not installed this returns false for all inputs.
"""
function is_darray(A)
    if !supports_distributed()
        return false
    end
    try
        import DistributedArrays
    catch
        return false
    end
    return A isa DistributedArrays.DArray
end

"""Gather a (possibly distributed) array on the `root` rank.

Behavior:
- If `DistributedArrays` is not installed then `A` is returned unchanged.
- If `A` is a `DArray` and the calling rank equals `root`, a full `Array` is
  returned containing the gathered data. On non-root ranks `nothing` is
  returned.

Arguments:
- `A`: array-like object (may be a `DArray`)
- `root`: integer rank that should receive the gathered array (default 0)
- `comm`: MPI communicator to use (default `MPI.COMM_WORLD`)
"""
function gather_on_root(A; root::Integer = 0, comm = MPI.COMM_WORLD)
    if !supports_distributed()
        return A
    end
    try
        import DistributedArrays
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
