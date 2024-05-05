var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [MUMPS]","category":"page"},{"location":"api/#MUMPS.MUMPS","page":"API","title":"MUMPS.MUMPS","text":"module MUMPS\n\nBoth low-level interface with MUMPS 5.6.2 parallel direct solver C-library as well as convenient wrappers for some common uses for MUMPS.\n\nThe central work is done by the Mumps struct, which mirrors the internal structure used in MUMPS. Manipulations can be done directly on this object and then passed to Mumps via the function invoke_mumps! This mode of operation gives the user complete control as described in the MUMPS manual, though it exposes unsafe operations, so beware.\n\nMore convenient are the use of the functions mumps_solve, mumps_factorize, mumps_det, mumps_schur_complement, and mumps_select_inv, which all have mutating counterparts (such as mumps_solve!). These can take matrices and right hand sides directly, so, for example, the equation A*x=y, solved in Base by x=A\\y or LinearAlbegra.ldiv!(x,A,y), can be solved with MUMPS as x=mumps_solve(A,y), or mumps_solve!(x,A,y).\n\nThe package also extends Base.det, Base.\\, LinearAlgebra.ldiv! and LinearAlgebra.inv to work with mumps objects.\n\nNote, unless working with the low-level interace, we discourage setting the JOB parameter manually, as this can lead to unsafe operation.\n\nThe goal is to give the advanced user low-level access to MUMPS, while simultaneously giving the ordinary user safe functions that grant access to most of what MUMPS has to offer.\n\n\n\n\n\n","category":"module"},{"location":"api/#MUMPS.default_cntl32","page":"API","title":"MUMPS.default_cntl32","text":"Default single precision real parameters\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.default_cntl64","page":"API","title":"MUMPS.default_cntl64","text":"Default double precision real parameters\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.default_icntl","page":"API","title":"MUMPS.default_icntl","text":"Default integer parameters.\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.mumps_definite","page":"API","title":"MUMPS.mumps_definite","text":"const mumps_definite\n\nConstant indicating that a symmetric definite matrix will be analyzed and factorized\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.mumps_symmetric","page":"API","title":"MUMPS.mumps_symmetric","text":"const mumps_symmetric\n\nConstant indicating that a general symmetric matrix will be analyzed and factorized\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.mumps_unsymmetric","page":"API","title":"MUMPS.mumps_unsymmetric","text":"const mumps_unsymmetric\n\nConstant indicating that a general unsymmetric matrix will be analyzed and factorized\n\n\n\n\n\n","category":"constant"},{"location":"api/#MUMPS.MUMPSException","page":"API","title":"MUMPS.MUMPSException","text":"mutable struct MUMPSException\n\nException type raised in case of error.\n\n\n\n\n\n","category":"type"},{"location":"api/#MUMPS.Mumps","page":"API","title":"MUMPS.Mumps","text":"Mirror of structre in [sdcz]mumps_c.h.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.finalize-Tuple{Mumps}","page":"API","title":"Base.finalize","text":"finalize(mumps)\n\nTerminate a Mumps instance.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.associate_matrix!-Tuple{Mumps, Integer, Vector, Vector, Vector}","page":"API","title":"MUMPS.associate_matrix!","text":"associate_matrix!(mumps, n, irow, jcol, vals)\n\nRegister the sparse matrix given in coordinate format with the Mumps object mumps. This function makes it possible to define the matrix on the host only. If the matrix is defined on all nodes, there is no need to use this function.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.associate_matrix!-Union{Tuple{TA}, Tuple{T}, Tuple{Mumps{T}, AbstractArray{TA}}} where {T, TA}","page":"API","title":"MUMPS.associate_matrix!","text":"associate_matrix!(mumps,A)\n\nRegister the square matrix A to a mumps object. It internally converts A to be consistent with the ICNTL[5] setting.\n\nIf needed, it tries to convert element type of A to be consistent with type of mumps, throwing a warning in this case.\n\nSee also: associate_rhs!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.associate_rhs!-Tuple{Mumps, AbstractMatrix}","page":"API","title":"MUMPS.associate_rhs!","text":"associate_rhs!(mumps, rhs)\n\nRegister a dense or sparse RHS matrix or vector rhs to a mumps object. It internally converts rhs to be consistent with the ICNTL[20] setting, and additionally allocates mumps.rhs according to the ICNTL[21] setting.\n\nIf needed, it tries to convert element type of rhs to be consistent with type of mumps, throwing a warning in this case.\n\nSee also: associate_matrix!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.default_icntl!-Tuple{Mumps}","page":"API","title":"MUMPS.default_icntl!","text":"default_icntl!(mumps)\n\nreset ICNTL to its default\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.display_cntl","page":"API","title":"MUMPS.display_cntl","text":"display_cntl(mumps)\n\nShow the complete CNTL real array of mumps, with descriptions\n\nSee also: set_cntl!\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.display_icntl","page":"API","title":"MUMPS.display_icntl","text":"display_icntl(mumps)\n\nShow the complete ICNTL integer array of mumps, with descriptions\n\nSee also: set_icntl!\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.factorize!-Tuple{Mumps, AbstractArray}","page":"API","title":"MUMPS.factorize!","text":"factorize!(mumps,A)\n\nCombined associate_matrix / factorize. Presume that A is available on all nodes.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.factorize!-Tuple{Mumps}","page":"API","title":"MUMPS.factorize!","text":"factorize!(mumps)\n\nFactorize the matrix registered with the Mumps instance. The matrix must have been previously registered with associate_matrix(). After the factorization, the determinant, if requested, is stored in mumps.det. The MUMPS error code is stored in mumps.err.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.finalize!-Tuple{Mumps}","page":"API","title":"MUMPS.finalize!","text":"finalize!(mumps)\n\nRelease the pointers contained in mumps\n\nSee also: initialize!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_icntl-Tuple{}","page":"API","title":"MUMPS.get_icntl","text":"get_icntl(;det=false,verbose=false,ooc=false,itref=0)\n\nObtain an array of integer control parameters.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_rhs!-Tuple{Any, Mumps}","page":"API","title":"MUMPS.get_rhs!","text":"get_rhs!(x,mumps)\n\nRetrieve right hand side from mumps, storing it in pre-allocated x\n\nSee also: get_rhs, get_sol!, get_sol\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_rhs-Union{Tuple{Mumps{T}}, Tuple{T}} where T","page":"API","title":"MUMPS.get_rhs","text":"get_rhs(mumps) -> y\n\nRetrieve right hand side from mumps\n\nSee also: get_rhs!, get_sol!, get_sol\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_schur_complement!-Tuple{Any, Mumps}","page":"API","title":"MUMPS.get_schur_complement!","text":"get_schur_complement!(S,mumps)\n\nRetrieve Schur complement matrix from mumps into pre-allocated S\n\nSee also: get_schur_complement, mumps_schur_complement!, mumps_schur_complement\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_schur_complement-Union{Tuple{Mumps{T}}, Tuple{T}} where T","page":"API","title":"MUMPS.get_schur_complement","text":"get_schur_complement(mumps) -> S\n\nRetrieve Schur complement matrix S from mumps\n\nSee also: get_schur_complement!, mumps_schur_complement!, mumps_schur_complement\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_sol!-Tuple{Union{SubArray, Array}, Mumps}","page":"API","title":"MUMPS.get_sol!","text":"get_sol!(x,mumps)\n\nRetrieve solution x from mumps into pre-allocated array.\n\nSee also: get_rhs!, get_rhs, get_sol\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_sol-Union{Tuple{Mumps{T}}, Tuple{T}} where T","page":"API","title":"MUMPS.get_sol","text":"get_sol(mumps) -> x\n\nRetrieve solution from mumps\n\nSee also: get_rhs!, get_rhs, get_sol!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.get_solution","page":"API","title":"MUMPS.get_solution","text":"get_solution(mumps) -> x\n\nRetrieve the solution of the system solved by solve(). This function makes it possible to ask MUMPS to assemble the final solution on the host only, and to retrieve it there.\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.initialize!-Tuple{Mumps}","page":"API","title":"MUMPS.initialize!","text":"initialize!(mumps)\n\nReinitialize mumps, regardless of its current state\n\nSee also: finalize\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.invoke_mumps!-Tuple{Mumps}","page":"API","title":"MUMPS.invoke_mumps!","text":"invoke_mumps!(mumps)\n\nCall the appropriate mumps C-library, passing to it the Mumps object mumps, but checking to make sure mumps has been initialized first, so that it's safe.\n\nThis is a low-level function, meaning that you have complete control over what operations are done, based on the MUMPS manual.\n\nBe warned, a direct call can crash Julia if mumps is not appropriately initialized.\n\nSee also: invoke_mumps_unsafe!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.invoke_mumps_unsafe!","page":"API","title":"MUMPS.invoke_mumps_unsafe!","text":"invoke_mumps_unsafe!(mumps)\n\nCall the appropriate mumps C-library, passing to it the Mumps object mumps.\n\nThis is a low-level function, meaning that you have complete control over what operations are done, based on the MUMPS manual.\n\nBe warned, a direct call can crash Julia if mumps is not appropriately initialized.\n\nSee also: invoke_mumps!\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.mumps_det!-Tuple{Mumps}","page":"API","title":"MUMPS.mumps_det!","text":"mumps_det!(mumps; discard=true)\n\nCompute determinant of A, which has been previously provided to mumps.\n\nDeterminant can be computed from mutated mumps by just det(mumps) [must have loaded LinearAlgebra].\n\nOptional keyward discard controls whether LU factors are discarded via ICNTL[31]. This is useful if you only care about the determinant and don't want to do any further computation with mumps. Use discard=2 to throw away only L.\n\nSee also: mumps_det\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.mumps_det-Tuple{Any}","page":"API","title":"MUMPS.mumps_det","text":"mumps_det(A) -> det\n\nCompute determinant of A.\n\nSee also: mumps_det!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.mumps_factorize!-Tuple{Mumps}","page":"API","title":"MUMPS.mumps_factorize!","text":"mumps_factorize!(mumps)\n\nLU factorize A previously provided to mump. LU stored in mumps, but not in a particularly accessible way. Useful for doing repeated solves downstream.\n\nSee also: mumps_factorize\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.mumps_factorize-Tuple{AbstractArray}","page":"API","title":"MUMPS.mumps_factorize","text":"mumps_factorize(A) -> mumps\n\nLU factorize A. LU stored in mumps, but not in a particularly accessible way. Useful for doing repeated solves downstream.\n\nSee also: mumps_factorize!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.mumps_schur_complement!","page":"API","title":"MUMPS.mumps_schur_complement!","text":"mumps_schur_complement!(mumps, schur_inds)\nmumps_schur_complement!(mumps, x)\n\nschur_inds is integer array of Schur indices. If x is sparse, Schur indices determined from populated rows of x\n\nSee also: mumps_schur_complement, get_schur_complement!, get_schur_complement\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.mumps_schur_complement-Tuple{AbstractArray, Any}","page":"API","title":"MUMPS.mumps_schur_complement","text":"mumps_schur_complement(A,schur_inds) -> S\nmumps_schur_complement(A,x) -> S\n\nschur_inds is integer array x is sparse, populated rows are Schur indices S is Schur complement matrix.\n\nSee also: mumps_schur_complement!\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.mumps_select_inv","page":"API","title":"MUMPS.mumps_select_inv","text":"mumps_select_inv(A,x) -> A⁻¹\nmumps_select_inv(A,I,J) -> A⁻¹\n\nCompute selected elements of A⁻¹. If two arguments are passed, the second must be sparse, and its sparsity pattern determines the entries of the inverse. If three arguments are passed, the integer arrays I and J specify which entries via i(k),j(k) = I[k],J[k].\n\nSee also: mumps_select_inv!\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.mumps_select_inv!","page":"API","title":"MUMPS.mumps_select_inv!","text":"mumps_select_inv!(x,mumps)\nmumps_select_inv!(x,A)\n\nCompute selected elements of A⁻¹ with same sparsity pattern as x, stored in x. If passed mumps, must have previously been provided with matrix A.\n\nSee also: mumps_select_inv, get_rhs!, get_rhs\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.mumps_solve","page":"API","title":"MUMPS.mumps_solve","text":"mumps_solve(A,y) -> x\nmumps_solve(mumps,y) -> x\nmumps_solve(mumps) -> x\n\nSolve A*x=y If mumps is given, must have previously been provided a matrix A. If only input is mumps must also have been provided y.\n\nSee also: mumps_solve!\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.mumps_solve!","page":"API","title":"MUMPS.mumps_solve!","text":"mumps_solve!(x,A,y; kwargs...)\nmumps_solve!(x,mumps)\nmumps_solve!(x,mumps,y)\n\nSolve A*x=y, saving result in pre-allocated x. If mumps is given, must have previously been provided a matrix A. If y is not given, mumps must have previously been provided y\n\nSee also: mumps_solve, get_sol!, get_sol\n\n\n\n\n\n","category":"function"},{"location":"api/#MUMPS.set_cntl!-Union{Tuple{TR}, Tuple{TC}, Tuple{Mumps{TC, TR}, Integer, AbstractFloat}} where {TC, TR}","page":"API","title":"MUMPS.set_cntl!","text":"set_cntl!(mumps,i,val; [displaylevel=1])\n\nSet the real/complex control parameters according to CNTL[i]=val\n\nSee also: display_cntl\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.set_icntl!-Tuple{Mumps, Integer, Integer}","page":"API","title":"MUMPS.set_icntl!","text":"set_icntl!(mumps,i,val; [displaylevel=1])\n\nSet the integer control parameters according to ICNTL[i]=val\n\nSee also: display_icntl\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.set_job!-Tuple{Mumps, Any}","page":"API","title":"MUMPS.set_job!","text":"set_job!(mumps,job)\n\nSet the phase to job. See MUMPS manual for options.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.set_save_dir!-Tuple{Any, String}","page":"API","title":"MUMPS.set_save_dir!","text":"set_save_dir!(mumps,dir)\n\nset name of directory in which to store out-of-core files.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.set_save_prefix!-Tuple{Any, String}","page":"API","title":"MUMPS.set_save_prefix!","text":"set_save_prefix!(mumps,dir)\n\nprefix for out-of-core files.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.set_schur_centralized_by_column!-Union{Tuple{T}, Tuple{Mumps{T}, AbstractArray{Int64}}} where T","page":"API","title":"MUMPS.set_schur_centralized_by_column!","text":"set_schur_centralized_by_column!(mumps,schur_inds)\n\nSet up Schur complement matrix calculation for the \"centralized by column\" method suggested in the MUMPS manual\n\nSee also: mumps_schur_complement!, mumps_schur_complement\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.solve!-Tuple{Mumps}","page":"API","title":"MUMPS.solve!","text":"solve!(mumps;transposed=false)\n\nSolve the system registered with the Mumps object mumps. The matrix and right-hand side(s) must have been previously registered with associate_matrix() and associate_rhs(). The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is stored internally and must be retrieved with get_solution().\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.solve-Tuple{Mumps, AbstractArray, AbstractArray}","page":"API","title":"MUMPS.solve","text":"solve(mumps, A, rhs; transposed=false)\n\nCombined analyze / factorize / solve. Presume that A and rhs are available on all nodes. The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is retrieved and returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.solve-Tuple{Mumps, AbstractArray}","page":"API","title":"MUMPS.solve","text":"solve(mumps, rhs; transposed=false)\n\nCombined associate_rhs / solve. Presume that rhs is available on all nodes. The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is retrieved and returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#MUMPS.solve-Union{Tuple{V}, Tuple{T}, Tuple{AbstractArray{T}, AbstractArray{V}}} where {T, V}","page":"API","title":"MUMPS.solve","text":"solve(A, rhs; sym=mumps_unsymmetric)\n\nCombined initialize / analyze / factorize / solve. Presume that A and rhs are available on all nodes. The optional keyword argument sym indicates the symmetry of A. The solution is retrieved and returned.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"#Home","page":"Home","title":"MUMPS.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MUMPS is a library for the solution of large linear systems using a factorization. Structure can be exploited, such as symmetry, or symmetry and definiteness. The factorization and solve phases can be performed in parallel via MPI by way of MPI.jl.","category":"page"},{"location":"#How-to-Install","page":"Home","title":"How to Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> ]\npkg> add MUMPS\npkg> test MUMPS","category":"page"},{"location":"#How-to-Use","page":"Home","title":"How to Use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main data type holding information on a factorization is Mumps. Remember to initialize MPI before attempting to create a Mumps object. A simple session is as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using MUMPS, MPI, SparseArrays, LinearAlgebra\nMPI.Init()\nA = sprand(10, 10, 0.2) + I\nrhs = rand(10)\nx = solve(A, rhs)\nnorm(x - A \\ rhs) / norm(x)\nMPI.Finalize()","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is possible to separate the initialization, the analysis/factorization, and the solve phases. It is also possible to access the information reported by MUMPS after the factorization and solve phases, and to modify this information (e.g., to perform iterative refinement).","category":"page"},{"location":"","page":"Home","title":"Home","text":"When creating an instance of a Mumps object explicitly, it is important to specify in advance what arithmetic should be used. Single and double precision real (Float32 and Float64) and complex (ComplexF32 and ComplexF64) arithmetics are supported.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For instance,","category":"page"},{"location":"","page":"Home","title":"Home","text":"MPI.Init()\nmumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64)\nA = sparse(rand(4,4))\nrhs = rand(4)\nassociate_matrix!(mumps, A)\nfactorize!(mumps)\nassociate_rhs!(mumps, rhs)\nsolve!(mumps)\nx = get_solution(mumps)\nfinalize(mumps)\nMPI.Finalize()","category":"page"},{"location":"","page":"Home","title":"Home","text":"Once the arithmetic of the Mumps instance has been specified, it cannot be changed. The module is flexible in that various data types may be used to define the matrix to be factorized and the right-hand side, and appropriate conversions will take place. Dense matrices may be used, and they will be converted to sparse format.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For intance,","category":"page"},{"location":"","page":"Home","title":"Home","text":"mumps = Mumps{ComplexF64}(mumps_unsymmetric, default_icntl, default_cntl64)\nA = rand(Int16, 4, 4)\nrhs = rand(Float32, 4)\nassociate_matrix!(mumps, A)  # A is converted to a sparse ComplexF64 matrix\nassociate_rhs!(mumps, rhs)   # rhs is converted to a Complex64 vector","category":"page"},{"location":"","page":"Home","title":"Home","text":"See test for more examples.","category":"page"},{"location":"#Constants-and-Methods-Exposed","page":"Home","title":"Constants and Methods Exposed","text":"","category":"section"},{"location":"#Constants","page":"Home","title":"Constants","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following convenience constants may be used when initializing a Mumps object:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Constant Meaning\nmumps_unsymmetric matrix is general unsymmetric (or symmetry is unknown)\nmumps_definite matrix is symmetric and (positive or negative) definite\nmumps_symmetric matrix is symmetric but indefinite (or definiteness is unknown)\ndefault_icntl array of default integer parameters\ndefault_cntl32 array of default real parameters in single precision\ndefault_cntl64 array of default real parameters in double precision","category":"page"},{"location":"","page":"Home","title":"Home","text":"See Section 5 of the MUMPS User's Manual for a description of the integer and real control arrays.","category":"page"},{"location":"#Methods","page":"Home","title":"Methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Mumps object is created using the default constructor, which must be supplied with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"the data type for the arithmetic to be used, as a type parameter, i.e., Mumps{Float64}(...) or Mumps{ComplexF64}(...)\nsym: one of the constants mumps_unsymmetric, mumps_definite or mumps_symmetric. Note that there is no support for Hermitian complex matrices in MUMPS. Therefore, we recommend to always use mumps_unsymmetric for complex data.\nicntl: an integer parameters array (see the MUMPS Users's Manual)\ncntl: a real parameters array (see the MUMPS Users's Manual)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The convenience function get_icntl() returns an array of integer parameters corresponding to certain commonly-used options. Its arguments are all optional:","category":"page"},{"location":"","page":"Home","title":"Home","text":"det: a boolean indicating whether the determinant should be computed (default: false)\nverbose: a boolean (default: false)\nooc: a boolean indicating whether factors should be stored out of core (default: false)\nitref: the number of iterative refinement steps (default: 0).","category":"page"},{"location":"","page":"Home","title":"Home","text":"A Mumps object is destroyed by calling the finalize() method. Because finalize still issues MPI commands, it is important to call finalize() before calling MPI.Finalize().","category":"page"},{"location":"","page":"Home","title":"Home","text":"Method Description\nfinalize Finalize a Mumps object. Must be done before calling MPI.Finalize()\nassociate_matrix! Register a matrix with the Mumps object. This function makes it possible to define the data on the host only.\nfactorize! Factorize the matrix registered with the Mumps object.\nassociate_rhs! Register right-hand sides with the Mumps object. This function makes it possible to define the data on the host only.\nsolve! Solve the linear system for the given right-hand side.\nget_solution Retrieve the solution from the Mumps object. This function makes it possible for the solution to be assembled on the host only.","category":"page"},{"location":"#Parallel-Execution","page":"Home","title":"Parallel Execution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MPI is controled by way of MPI.jl. Look for the lines that say NUMBER OF WORKING PROCESSES in the output of","category":"page"},{"location":"","page":"Home","title":"Home","text":"mpirun -np 4 julia examples/mumps_mpi.jl","category":"page"},{"location":"#Custom-Installation","page":"Home","title":"Custom Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Note: MUMPS is already precompiled with Yggdrasil for all platforms except Windows.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To use your custom MUMPS, set the environment variable JULIA_MUMPS_LIBRARY_PATH to point to the shared library before using MUMPS. Note that the same version of MUMPS as used by the MUMPS_jll artifact is needed.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, macOS users may install precompiled MUMPS binaries from the Homebrew tap dpo/mumps-jl as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"brew tap dpo/mumps-jl\nbrew install mpich-mumps\nexport JULIA_MUMPS_LIBRARY_PATH=$(brew --prefix)/opt/mpich-mumps/lib","category":"page"},{"location":"","page":"Home","title":"Home","text":"Apple Silicon users should remember to use arch x86_64 brew to refer to Intel binaries run through Rosetta, as we do not (yet) ship Silicon binaries of MUMPS via Homebrew.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The JULIA_MUMPS_LIBRARY_PATH environment variable may be set permanently in the shell's startup file, or in $HOME/.julia/config/startup.jl.","category":"page"}]
}
