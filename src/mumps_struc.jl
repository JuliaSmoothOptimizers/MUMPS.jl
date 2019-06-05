# this file mirros the relevant content of the "[sdcz]mumps_c.h" file of MUMPS 5.2.0
# `gc_haven`, contains Julia references to protect the pointers passed to C
# from gargage collection.
export Mumps, MUMPSException

"""
    mutable struct MUMPSException
Exception type raised in case of error.
"""
mutable struct MUMPSException <: Exception
  msg :: AbstractString
end
const MUMPSValueDataType = Union{Float32, Float64, ComplexF32, ComplexF64}

const MUMPS_VERSION = "5.2.0"
const MUMPS_VERSION_MAX_LEN = 30
const DEFAULT_FORTRAN_COMMUNICATOR = -987654

# mirror of structre in [sdcz]mumps_c.h
mutable struct Mumps{TC,TR}
    sym::MUMPS_INT # MANDATORY 0 for unsymmetric, 1 for symmetric and posdef, 2 for general symmetric. All others treated as 0
    par::MUMPS_INT # MANDATORY 0 host not involved in parallel factorization and solve, 1 host is involved
    job::MUMPS_INT # MANDATORY -1 initializes package, must come first, -2 terminates, 1 analysis, 2 factorization, 3 solve, 4=1&2, 5=2&3, 6=1&2&3
    comm_fortran::MUMPS_INT # MANDATORY valid MPI communicator
    icntl::NTuple{60,MUMPS_INT}
    keep::NTuple{500,MUMPS_INT}
    cntl::NTuple{15,TR}
    dkeep::NTuple{230,TR}
    keep8::NTuple{150,MUMPS_INT8}
    n::MUMPS_INT

    nz_alloc::MUMPS_INT

    nz::MUMPS_INT
    nnz::MUMPS_INT8
    irn::Ptr{MUMPS_INT}
    jcn::Ptr{MUMPS_INT}
    a::Ptr{TC}

    nz_loc::MUMPS_INT
    nnz_loc::MUMPS_INT8
    irn_loc::Ptr{MUMPS_INT}
    jcn_loc::Ptr{MUMPS_INT}
    a_loc::Ptr{TC}

    nelt::MUMPS_INT
    eltptr::Ptr{MUMPS_INT}
    eltvar::Ptr{MUMPS_INT}
    a_elt::Ptr{TC}

    perm_in::Ptr{MUMPS_INT}

    sym_perm::Ptr{MUMPS_INT}
    uns_perm::Ptr{MUMPS_INT}

    colsca::Ptr{TR}
    rowsca::Ptr{TR}

    colsca_from_mumps::MUMPS_INT
    rowsca_from_mumps::MUMPS_INT

    rhs::Ptr{TC}
    redrhs::Ptr{TC}
    rhs_sparse::Ptr{TC}
    sol_loc::Ptr{TC}
    rhs_loc::Ptr{TC}

    irhs_sparse::Ptr{MUMPS_INT}
    irhs_ptr::Ptr{MUMPS_INT}
    isol_loc::Ptr{MUMPS_INT}
    irhs_loc::Ptr{MUMPS_INT}

    nrhs::MUMPS_INT
    lrhs::MUMPS_INT
    lredrhs::MUMPS_INT
    nz_rhs::MUMPS_INT
    lsol_loc::MUMPS_INT
    nloc_rhs::MUMPS_INT
    lrhs_loc::MUMPS_INT

    schur_mloc::MUMPS_INT
    schur_nloc::MUMPS_INT
    schur_lld::MUMPS_INT

    mblock::MUMPS_INT
    nblock::MUMPS_INT
    nprow::MUMPS_INT
    npcol::MUMPS_INT

    info::NTuple{80,MUMPS_INT}
    infog::NTuple{80,MUMPS_INT}

    rinfo::NTuple{40,TR}
    rinfog::NTuple{40,TR}

    deficiency::MUMPS_INT
    pivnul_list::Ptr{MUMPS_INT}
    mapping::Ptr{MUMPS_INT}

    size_schur::MUMPS_INT
    listvar_schur::Ptr{MUMPS_INT}
    schur::Ptr{TC}

    instance_number ::MUMPS_INT
    wk_user     ::Ptr{TC}

    version_number ::NTuple{MUMPS_VERSION_MAX_LEN+1+1,Cchar}
    ooc_tmpdir     ::NTuple{256,Cchar}
    ooc_prefix      ::NTuple{64,Cchar}
    write_problem  ::NTuple{256,Cchar}
    lwk_user    ::MUMPS_INT
    save_dir    ::NTuple{256,Cchar}
    save_prefix ::NTuple{256,Cchar}

    metis_options::NTuple{40,MUMPS_INT}

    det::TC
    err::Int

    _gc_haven::Array{Ref,1}
    _finalized::Bool

    function Mumps{T}(sym::Integer,par::Integer,comm::Integer) where T<:MUMPSValueDataType
        !MPI.Initialized() ? throw(MUMPSException("Initialize MPI first")) : nothing
        mumps = new{T,real(T)}(sym,par,-1,comm)
        invoke_mumps_unsafe!(mumps)
        mumps._gc_haven = Array{Ref,1}(undef,0)
        mumps._finalized = false
        finalizer(finalize!,mumps)
        return mumps
    end
end

Mumps{T}(sym::Integer,par::Integer=1) where T<:MUMPSValueDataType = Mumps{T}(sym,par,DEFAULT_FORTRAN_COMMUNICATOR)
function Mumps{T}(sym::Integer,par::Integer=1) where T<:Number
    if promote_type(T,Float32)<:MUMPSValueDataType
        Mumps{promote_type(T,Float32)}(sym,par)
    elseif promote_type(T,Float64)<:MUMPSValueDataType
        Mumps{promote_type(T,Float32)}(sym,par)
    elseif promote_type(T,ComplexF32)<:MUMPSValueDataType
        Mumps{promote_type(T,Float32)}(sym,par)
    elseif promote_type(T,ComplexF64)<:MUMPSValueDataType
        Mumps{promote_type(T,Float32)}(sym,par)
    else
        throw(MUMPSException("cannot promote type $T to a MUMPS-compatible type"))
    end
end

function Mumps{T}(
            sym::Integer,
            icntl::Array{TI,1},
            cntl::Array{V,1}; par=1
            ) where {TI<:Integer, T<:MUMPSValueDataType,V<:AbstractFloat}

    # Set default pivot threshold if required.
    if cntl[1] == -1
          cntl[1] = (sym == mumps_definite) ? 0.0 : 0.01
    end
    mumps = Mumps{T}(sym, par)
    for i ∈ eachindex(icntl)
        set_icntl!(mumps,i,icntl[i]; displaylevel=0)
    end
    for i ∈ eachindex(cntl)
        set_cntl!(mumps,i,cntl[i]; displaylevel=0)
    end
    return mumps
end

# this is necessary to make the call "real(TV)" work in the inner constructor for Mumps above
function Base.real(T::TypeVar)
    if T<:Number
        return real(T)
    else
        throw(DomainError("real not defined for type $T"))
    end
end
