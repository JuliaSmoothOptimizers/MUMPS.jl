"""
    module MUMPS

Both low-level interface with MUMPS 5.2.0 parallel direct solver C-library
as well as convenient wrappers for some common uses for MUMPS.

The central work is done by the `Mumps` struct, which mirrors the
internal structure used in MUMPS. Manipulations can be done directly
on this object and then passed to Mumps via the function [`invoke_mumps!`](@ref)
This mode of operation gives the user complete control as described
in the MUMPS manual, though it exposes unsafe operations, so beware.

More convenient are the use of the functions [`mumps_solve`](@ref), [`mumps_factorize`](@ref),
[`mumps_det`](@ref), [`mumps_schur`](@ref), and [`mumps_select_inv`](@ref), which all have
mutating counterparts (such as [`mumps_solve!`](@ref)). These can take matrices
and right hand sides directly, so, for example, the equation `A*x=y`, solved
in Base by `x=A\\y` or `LinearAlbegra.ldiv!(x,A,y)`, can be solved in MUMPS3
as `x=mumps_solve(A,y)`, or `mumps_solve!(x,A,y)`.

The package also extends Base.det, Base.\\, LinearAlgebra.ldiv! and LinearAlgebra.inv to
work with mumps objects.

Note, unless working with the low-level interace, I discourage setting the `JOB`
parameter manually, as this can lead to unsafe operation.

The goal is to give the advanced user low-level access to MUMPS, while simultaneously
giving the ordinary user safe functions that grant access to most of what
MUMPS has to offer.

(https://github.com/wrs28/MUMPS3.jl).
"""
# module MUMPS

using MPI,
Libdl,
LinearAlgebra,
SparseArrays

function __init__()
  if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
  else
    error("MUMPS library not properly installed. Please run Pkg.build(\"MUMPS\")")
  end
  global LIB_S = dlopen(libsmumps)
  global LIB_D = dlopen(libdmumps)
  global LIB_C = dlopen(libcmumps)
  global LIB_Z = dlopen(libzmumps)
end

include("mumps_types.jl")
include("mumps_struc.jl")
include("interface.jl")
include("convenience.jl")
include("icntl_alibis.jl")
include("printing.jl")

include("exported_methods.jl")

# end
