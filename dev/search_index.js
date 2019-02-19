var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Home-1",
    "page": "Home",
    "title": "MUMPS.jl documentation",
    "category": "section",
    "text": "MUMPS is a library for the solution of large linear systems using a factorization. Structure can be exploited, such as symmetry, or symmetry and definiteness. The factorization and solve phases can be performed in parallel."
},

{
    "location": "#How-to-Install-1",
    "page": "Home",
    "title": "How to Install",
    "category": "section",
    "text": ""
},

{
    "location": "#Prerequisites-1",
    "page": "Home",
    "title": "Prerequisites",
    "category": "section",
    "text": "Because BinDeps is essentially broken, you must install MUMPS outside of Julia. On macOS, we recommend using Homebrew. On Linux, we recommend using Linuxbrew. In both cases, the commands are the same:$ brew tap brewsci/num\n$ brew install brewsci-mumps  # use brew options brewsci-mumps for build optionsNote: apt-get install libmumps-dev installs a version of OpenMPI that is too old for MPI.jl. See the Troubleshooting section below.All examples above install OpenMPI. If you wish to use MPICH, you will have to build MUMPS by hand."
},

{
    "location": "#Building-MUMPS.jl-1",
    "page": "Home",
    "title": "Building MUMPS.jl",
    "category": "section",
    "text": "If MUMPS and SCALAPACK are not in standard locations where BinDeps will find them, you can help by setting the environment variables MUMPS_PREFIX and SCALAPACK_PREFIX.The Homebrew and Linuxbrew methods above install MUMPS and SCALAPACK in nonstandard locations. You can definejulia> ENV[\"MUMPS_PREFIX\"] = \"/usr/local/opt/brewsci-mumps\"\njulia> ENV[\"SCALAPACK_PREFIX\"] = \"/usr/local/opt/brewsci-scalapack\"on macOS, and something of the formjulia> ENV[\"MUMPS_PREFIX\"] = \"/home/linuxbrew/.linuxbrew/opt/brewsci-mumps\"\njulia> ENV[\"SCALAPACK_PREFIX\"] = \"/home/linuxbrew/.linuxbrew/opt/brewsci-scalapack\"on Linux.At the Julia prompt, typejulia> using Pkg\njulia> Pkg.clone(\"https://github.com/JuliaSmoothOptimizers/MUMPS.jl.git\")\njulia> Pkg.build(\"MUMPS\")\njulia> Pkg.test(\"MUMPS\")"
},

{
    "location": "#Troubleshooting-1",
    "page": "Home",
    "title": "Troubleshooting",
    "category": "section",
    "text": "On macOS or Linux, if you see the error message[ 11%] Building Fortran object CMakeFiles/gen_constants.dir/gen_constants.f90.o\n│ /home/ubuntu/.julia/packages/MPI/U5ujD/deps/gen_constants.f90:43:43:\n│\n│    call output(\"MPI_NO_OP       \", MPI_NO_OP)\n│                                            1\n│ Error: Symbol ‘mpi_no_op’ at (1) has no IMPLICIT typeyour OpenMPI library is too old.If you are running macOS and see error messages of the formPMIx has detected a temporary directory name that results in a path that is too long for the Unix domain socket:\n\n  Temp dir:\n  /var/folders/rq/p5nq9tv17p5drlk49755jjz80000gn/T/openmpi-sessions-501@your_computer_name_0/44473\n\nTry setting your TMPDIR environmental variable to point to something shorter in lengthsimply exit Julia and set the environment variable TMPDIR to, e.g., \\tmp:$ export TMPDIR=/tmpThe issue has to do with OpenMPI and is documented in their faq."
},

{
    "location": "#How-to-Use-1",
    "page": "Home",
    "title": "How to Use",
    "category": "section",
    "text": "The main data type holding information on a factorization is Mumps. Remember to initialize MPI before attempting to create a Mumps object. A simple session is as follows:julia> using MUMPS\njulia> using MPI\njulia> MPI.Init()\njulia> A = sprand(10, 10, .2) + speye(10); rhs = rand(10)\njulia> x = solve(A, rhs)  # Mumps object is created and destroyed\njulia> norm(x - A \\ rhs) / norm(x)\n2.640677159735313e-16\njulia> MPI.Finalize()     # if you\'re finishedIt is possible to separate the initialization, the analysis/factorization, and the solve phases. It is also possible to access the information reported by MUMPS after the factorization and solve phases, and to modify this information (e.g., to perform iterative refinement).When creating an instance of a Mumps object explicitly, it is important to specify in advance what arithmetic should be used. Single and double precision real (Float32 and Float64) and complex (Complex64 and Complex128) arithmetics are supported.For instance,julia> MPI.Init()\njulia> mumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64)  # Real, general unsymmetric\njulia> A = sparse(rand(4,4)); rhs = rand(4)       # Happens on all cores\njulia> associate_matrix!(mumps, A)\njulia> factorize!(mumps)\njulia> associate_rhs!(mumps, rhs)\njulia> solve!(mumps)\njulia> x = get_solution(mumps)\njulia> finalize(mumps)\njulia> MPI.Finalize()Once the arithmetic of the Mumps instance has been specified, it cannot be changed. The module is flexible in that various data types may be used to define the matrix to be factorized and the right-hand side, and appropriate conversions will take place. Dense matrices may be used, and they will be converted to sparse format.For intance,julia> mumps = Mumps{Complex128}(mumps_unsymmetric, default_icntl, default_cntl64)\njulia> A = rand(Int16, 4, 4); rhs = rand(Float32, 4)\njulia> associate_matrix!(mumps, A)  # A is converted to a sparse Complex128 matrix\njulia> associate_rhs!(mumps, rhs)   # rhs is converted to a Complex128 arraySee test for more examples."
},

{
    "location": "#Constants-and-Methods-Exposed-1",
    "page": "Home",
    "title": "Constants and Methods Exposed",
    "category": "section",
    "text": ""
},

{
    "location": "#Constants-1",
    "page": "Home",
    "title": "Constants",
    "category": "section",
    "text": "The following convenience constants may be used when initializing a Mumps object:Constant Meaning\nmumps_unsymmetric matrix is general unsymmetric (or symmetry is unknown)\nmumps_definite matrix is symmetric and (positive or negative) definite\nmumps_symmetric matrix is symmetric but indefinite (or definiteness is unknown)\ndefault_icntl array of default integer parameters\ndefault_cntl32 array of default real parameters in single precision\ndefault_cntl64 array of default real parameters in double precisionSee Sections 5.1 and 5.2 of the MUMPS User\'s Manual for a description of the integer and real control arrays."
},

{
    "location": "#Methods-1",
    "page": "Home",
    "title": "Methods",
    "category": "section",
    "text": "A Mumps object is created using the default constructor, which must be supplied with:the data type for the arithmetic to be used, as a type parameter, i.e., Mumps{Float64}(...) or Mumps{Complex128}(...)\nsym: one of the constants mumps_unsymmetric, mumps_definite or mumps_symmetric. Note that there is no support for Hermitian complex matrices in MUMPS. Therefore, we recommend to always use mumps_unsymmetric for complex data.\nicntl: an integer parameters array (see the MUMPS Users\'s Manual)\ncntl: a real parameters array (see the MUMPS Users\'s Manual)The convenience function get_icntl() returns an array of integer parameters corresponding to certain commonly-used options. Its arguments are all optional:det: a boolean indicating whether the determinant should be computed (default: false)\nverbose: a boolean (default: false)\nooc: a boolean indicating whether factors should be stored out of core (default: false)\nitref: the number of iterative refinement steps (default: 0).A Mumps object is destroyed by calling the finalize() method. Because finalize still issues MPI commands, it is important to call finalize() before calling MPI.Finalize().Method Description\nfinalize Finalize a Mumps object. Must be done before calling MPI.Finalize()\nassociate_matrix! Register a matrix with the Mumps object. This function makes it possible to define the data on the host only.\nfactorize! Factorize the matrix registered with the Mumps object.\nassociate_rhs! Register right-hand sides with the Mumps object. This function makes it possible to define the data on the host only.\nsolve! Solve the linear system for the given right-hand side.\nget_solution Retrieve the solution from the Mumps object. This function makes it possible for the solution to be assembled on the host only."
},

{
    "location": "#Parallel-Execution-1",
    "page": "Home",
    "title": "Parallel Execution",
    "category": "section",
    "text": "MPI is controled by way of MPI.jl. Look for the lines that say NUMBER OF WORKING PROCESSES in the output ofmpirun -np 4 julia examples/mumps_mpi.jl"
},

{
    "location": "#To-Do-(Pull-Requests-Welcome!)-1",
    "page": "Home",
    "title": "To Do (Pull Requests Welcome!)",
    "category": "section",
    "text": "[X] Support double precision complex arithmetic (in 99c23fe)\n[X] Support single precision real and complex arithmetic (in 654814a)\n[ ] Support distributed matrices / vectors\n[ ] User-selected permutation\n[X] Out-of-core option (in 73e829b)\n[X] Determinant (in 73e829b)\n[ ] Compute entries of the inverse\n[X] Control iterative refinement (in 73e829b)\n[ ] Obtain a Schur complement\n[ ] Solve with sparse right-hand sides\n[ ] Sequential, version with no MPI requirement"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Pages = [\"api.md\"]"
},

{
    "location": "api/#MUMPS.default_icntl",
    "page": "API",
    "title": "MUMPS.default_icntl",
    "category": "constant",
    "text": "Default integer parameters.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.default_cntl32",
    "page": "API",
    "title": "MUMPS.default_cntl32",
    "category": "constant",
    "text": "Default single precision real parameters\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.default_cntl64",
    "page": "API",
    "title": "MUMPS.default_cntl64",
    "category": "constant",
    "text": "Default double precision real parameters\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.mumps_unsymmetric",
    "page": "API",
    "title": "MUMPS.mumps_unsymmetric",
    "category": "constant",
    "text": "Constant indicating that a general unsymmetric matrix will be analyzed and factorized\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.mumps_definite",
    "page": "API",
    "title": "MUMPS.mumps_definite",
    "category": "constant",
    "text": "Constant indicating that a symmetric definite matrix will be analyzed and factorized\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.mumps_symmetric",
    "page": "API",
    "title": "MUMPS.mumps_symmetric",
    "category": "constant",
    "text": "Constant indicating that a general symmetric matrix will be analyzed and factorized\n\n\n\n\n\n"
},

{
    "location": "api/#Constants-1",
    "page": "API",
    "title": "Constants",
    "category": "section",
    "text": "default_icntl\ndefault_cntl32\ndefault_cntl64\nmumps_unsymmetric\nmumps_definite\nmumps_symmetric"
},

{
    "location": "api/#MUMPS.MUMPSException",
    "page": "API",
    "title": "MUMPS.MUMPSException",
    "category": "type",
    "text": "Exception type raised in case of error.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.Mumps",
    "page": "API",
    "title": "MUMPS.Mumps",
    "category": "type",
    "text": "Abstract type representing a factorization with MUMPS. All constructor arguments are optional. By default a general unsymmetric matrix will be analyzed/factorized with default integer and real parameters\n\n\n\n\n\n"
},

{
    "location": "api/#Types-1",
    "page": "API",
    "title": "Types",
    "category": "section",
    "text": "MUMPSException\nMumps"
},

{
    "location": "api/#MUMPS.get_icntl",
    "page": "API",
    "title": "MUMPS.get_icntl",
    "category": "function",
    "text": "Obtain an array of integer control parameters.\n\n\n\n\n\n"
},

{
    "location": "api/#Base.finalize",
    "page": "API",
    "title": "Base.finalize",
    "category": "function",
    "text": "Terminate a Mumps instance.\n\n\n\n\n\nTerminate a Mumps instance.\n\n\n\n\n\nTerminate a Mumps instance.\n\n\n\n\n\nTerminate a Mumps instance.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.factorize!",
    "page": "API",
    "title": "MUMPS.factorize!",
    "category": "function",
    "text": "Factorize the matrix registered with the Mumps instance. The matrix must have been previously registered with associate_matrix(). After the factorization, the determinant, if requested, is stored in mumps.det. The MUMPS error code is stored in mumps.err. \n\n\n\n\n\nFactorize the matrix registered with the Mumps instance. The matrix must have been previously registered with associate_matrix(). After the factorization, the determinant, if requested, is stored in mumps.det. The MUMPS error code is stored in mumps.err. \n\n\n\n\n\nFactorize the matrix registered with the Mumps instance. The matrix must have been previously registered with associate_matrix(). After the factorization, the determinant, if requested, is stored in mumps.det. The MUMPS error code is stored in mumps.err. \n\n\n\n\n\nFactorize the matrix registered with the Mumps instance. The matrix must have been previously registered with associate_matrix(). After the factorization, the determinant, if requested, is stored in mumps.det. The MUMPS error code is stored in mumps.err. \n\n\n\n\n\nCombined associate_matrix / factorize. Presume that A is available on all nodes.\n\n\n\n\n\nCombined associate_matrix / factorize. Presume that A is available on all nodes.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.solve!",
    "page": "API",
    "title": "MUMPS.solve!",
    "category": "function",
    "text": "Solve the system registered with the Mumps object mumps. The matrix and right-hand side(s) must have been previously registered with associate_matrix() and associate_rhs(). The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is stored internally and must be retrieved with get_solution().\n\n\n\n\n\nSolve the system registered with the Mumps object mumps. The matrix and right-hand side(s) must have been previously registered with associate_matrix() and associate_rhs(). The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is stored internally and must be retrieved with get_solution().\n\n\n\n\n\nSolve the system registered with the Mumps object mumps. The matrix and right-hand side(s) must have been previously registered with associate_matrix() and associate_rhs(). The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is stored internally and must be retrieved with get_solution().\n\n\n\n\n\nSolve the system registered with the Mumps object mumps. The matrix and right-hand side(s) must have been previously registered with associate_matrix() and associate_rhs(). The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is stored internally and must be retrieved with get_solution().\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.solve",
    "page": "API",
    "title": "MUMPS.solve",
    "category": "function",
    "text": "Combined associate_rhs / solve. Presume that rhs is available on all nodes. The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is retrieved and returned.\n\n\n\n\n\nCombined analyze / factorize / solve. Presume that A and rhs are available on all nodes. The optional keyword argument transposed indicates whether the user wants to solve the forward or transposed system. The solution is retrieved and returned.\n\n\n\n\n\nCombined initialize / analyze / factorize / solve. Presume that A and rhs are available on all nodes. The optional keyword argument sym indicates the symmetry of A. The solution is retrieved and returned.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.associate_matrix!",
    "page": "API",
    "title": "MUMPS.associate_matrix!",
    "category": "function",
    "text": "Register the matrix A with the Mumps object mumps. This function makes it possible to define the matrix on the host only. If the matrix is defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the matrix A with the Mumps object mumps. This function makes it possible to define the matrix on the host only. If the matrix is defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the matrix A with the Mumps object mumps. This function makes it possible to define the matrix on the host only. If the matrix is defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the matrix A with the Mumps object mumps. This function makes it possible to define the matrix on the host only. If the matrix is defined on all nodes, there is no need to use this function.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.associate_rhs!",
    "page": "API",
    "title": "MUMPS.associate_rhs!",
    "category": "function",
    "text": "Register the right-hand side(s) rhs with the Mumps object mumps. This function makes it possible to define the right- -hand side(s) on the host only. If the right-hand side(s) are defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the right-hand side(s) rhs with the Mumps object mumps. This function makes it possible to define the right- -hand side(s) on the host only. If the right-hand side(s) are defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the right-hand side(s) rhs with the Mumps object mumps. This function makes it possible to define the right- -hand side(s) on the host only. If the right-hand side(s) are defined on all nodes, there is no need to use this function.\n\n\n\n\n\nRegister the right-hand side(s) rhs with the Mumps object mumps. This function makes it possible to define the right- -hand side(s) on the host only. If the right-hand side(s) are defined on all nodes, there is no need to use this function.\n\n\n\n\n\n"
},

{
    "location": "api/#MUMPS.get_solution",
    "page": "API",
    "title": "MUMPS.get_solution",
    "category": "function",
    "text": "Retrieve the solution of the system solved by solve(). This function makes it possible to ask MUMPS to assemble the final solution on the host only, and to retrieve it there.\n\n\n\n\n\nRetrieve the solution of the system solved by solve(). This function makes it possible to ask MUMPS to assemble the final solution on the host only, and to retrieve it there.\n\n\n\n\n\nRetrieve the solution of the system solved by solve(). This function makes it possible to ask MUMPS to assemble the final solution on the host only, and to retrieve it there.\n\n\n\n\n\nRetrieve the solution of the system solved by solve(). This function makes it possible to ask MUMPS to assemble the final solution on the host only, and to retrieve it there.\n\n\n\n\n\n"
},

{
    "location": "api/#Utilities-1",
    "page": "API",
    "title": "Utilities",
    "category": "section",
    "text": "get_icntl\nfinalize\nfactorize!\nsolve!\nsolve\nassociate_matrix!\nassociate_rhs!\nget_solution"
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": ""
},

]}
