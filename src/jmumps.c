#include <stdio.h>
#include <stdlib.h>
#include "jmumps.h"

#define JOB_INIT                   -1
#define JOB_END                    -2
#define JOB_ANALYZE                 1
#define JOB_FACTORIZE               2
#define JOB_SOLVE                   3
#define JOB_ANALYZE_FACTORIZE       4
#define JOB_FACTORIZE_SOLVE         5
#define JOB_ANALYZE_FACTORIZE_SOLVE 6

#define USE_COMM_WORLD -987654

// Initialize MPI.
int mumps_initialize_mpi(void) {

  int    argc = 1;
  char** argv = NULL;
  int    taskid, np, ierr;

  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf("Error initializing MPI. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, ierr);
    return ierr;
  }

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  if (taskid == 0) printf("MPI initialized with %d processes\n", np);
  return ierr;
}

// Finalize MPI.
int mumps_finalize_mpi(void) {

  int ierr;

  ierr = MPI_Finalize();
  return ierr;
}

// Initialize
// Initialize MUMPS data structure, store default parameters.
void* mumps_initialize(int sym, int* icntl, double* cntl) {

  DMUMPS_STRUC_C* pmumps;
  int mpi_initialized, ierr;
  int i;

  MPI_Initialized(&mpi_initialized);
  if (!mpi_initialized) {
    printf("Please initialize MPI first\n");
    return NULL;
  }

#ifdef JMUMPS_DEBUG
  int taskid, np;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  printf("MPI process %d out of %d\n", np);
#endif

  // Initialize MUMPS.
  mumps_alloc(&pmumps);
  if (pmumps == NULL) {
    ierr = MPI_Finalize();
    return NULL;
  }

  pmumps->sym = sym;
  pmumps->job = JOB_INIT;
  pmumps->par = 1;
  pmumps->comm_fortran = USE_COMM_WORLD;

  pmumps->n   = -1;
  pmumps->nz  = -1;

  dmumps_c(pmumps);

  // Fill in default parameter values.
  for (i = 0; i < 40; i++)
    (pmumps->icntl)[i] = icntl[i];

  for (i = 0; i < 5; i++)
    (pmumps->cntl)[i] = cntl[i];

  return (void*)pmumps;
}

// Factorize
// Factorize input matrix.
void mumps_factorize(void *jmumps, int n, int nz,
                     double* vals, int* irow, int* jcol) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int taskid, ierr;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  // Analyze and factorize.
  if (taskid == 0) {

#ifdef JMUMPS_DEBUG
    printf("MUMPS struct initialized at %p\n", pmumps);
#endif

    pmumps->n  = n;
    pmumps->nz = nz;

    // irow/jcol are 1-based in Julia, as in Fortran.
    pmumps->irn = irow;
    pmumps->jcn = jcol;
    pmumps->a   = vals;

#ifdef JMUMPS_DEBUG
    printf("MUMPS factorized a matrix of size %d with %d nonzeros\n",
           pmumps->n, pmumps->nz);
#endif
  }

  pmumps->job = JOB_ANALYZE_FACTORIZE;
  dmumps_c(pmumps);
}

// Solve
// Solve a linear system using the factorization computed previously.
void mumps_solve(void* jmumps, int nrhs, double* rhs) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int taskid, ierr;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  if (taskid == 0) {

#ifdef JMUMPS_DEBUG
    printf("Solving linear system at %p\n", pmumps);
    printf("  nrhs=%d\n", nrhs);
    printf("     n=%d\n", pmumps->n);
    printf("    nz=%d\n", pmumps->nz);
    printf("   irn=%p\n", pmumps->irn);
    printf("   jcn=%p\n", pmumps->jcn);
    printf("     a=%p\n", pmumps->a);
    int k;
    for (k = 0; k < pmumps->n; k++) {
      printf("  rhs[%d]=%8.1e\n", k, rhs[k]);
    }
    for (k = 0; k < pmumps->nz; k++)
      printf("  irow[%d]=%d  jcol[%d]=%d  a[%d]=%8.1e\n",
             k, (pmumps->irn)[k], k, (pmumps->jcn)[k], k, (pmumps->a)[k]);
#endif

    pmumps->nrhs = nrhs;
    pmumps->lrhs = pmumps->n;
    pmumps->rhs  = rhs;  // Will be overwritten with the solution.
  }

  pmumps->job  = JOB_SOLVE;
  dmumps_c(pmumps);
  return;
}

// Get info
// Obtain analysis / factorization / solve integer and real info vectors.
void mumps_get_info(void* jmumps, int* infog, double* rinfog) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int taskid, ierr;
  int i;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  if (taskid == 0) {
    for (i = 0; i < 40; i++)
      infog[i] = (pmumps->infog)[i];

    for (i = 0; i < 20; i++)
      rinfog[i] = (pmumps->rinfog)[i];
  }
  return;
}

// Finalize
// Destroy MUMPS data structure.
void mumps_finalize(void* jmumps) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int taskid = 0, mpi_finalized = 0;

  if (pmumps == NULL) return;

  MPI_Finalized(&mpi_finalized);

  if (!mpi_finalized)
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

#ifdef JMUMPS_DEBUG
  if (taskid == 0)
    printf("Terminating MUMPS struct at %p\n", pmumps);
#endif

  pmumps->job = JOB_END;
  dmumps_c(pmumps);
  mumps_free(&pmumps);

#ifdef JMUMPS_DEBUG
  if (taskid == 0)
    printf("MUMPS instance terminated\n");
#endif

  return;
}


// Helper functions, strongly inspired by the MUMPS MATLAB interface.

void mumps_alloc(DMUMPS_STRUC_C** mumps){

  (*mumps) = malloc(sizeof(DMUMPS_STRUC_C));
  if (*mumps == NULL) return;

  (*mumps)->irn           = NULL;
  (*mumps)->jcn           = NULL;
  (*mumps)->a             = NULL;
  (*mumps)->irn_loc       = NULL;
  (*mumps)->jcn_loc       = NULL;
  (*mumps)->a_loc         = NULL;
  (*mumps)->eltptr        = NULL;
  (*mumps)->eltvar        = NULL;
  (*mumps)->a_elt         = NULL;
  (*mumps)->perm_in       = NULL;
  (*mumps)->colsca        = NULL;
  (*mumps)->rowsca        = NULL;
  (*mumps)->rhs           = NULL;
  (*mumps)->redrhs        = NULL;
  (*mumps)->rhs_sparse    = NULL;
  (*mumps)->irhs_sparse   = NULL;
  (*mumps)->irhs_ptr      = NULL;
  (*mumps)->pivnul_list   = NULL;
  (*mumps)->listvar_schur = NULL;
  (*mumps)->schur         = NULL;
  (*mumps)->sym_perm      = NULL;
  (*mumps)->uns_perm      = NULL;
}


#define Free(ptr) if (ptr) {free(ptr); ptr=NULL;}

void mumps_free(DMUMPS_STRUC_C** mumps){
  if (*mumps != NULL){
    // Free( (*mumps)->irn );
    // Free( (*mumps)->jcn  );
    // Free( (*mumps)->a );
    Free( (*mumps)->irn_loc );
    Free( (*mumps)->jcn_loc );
    Free( (*mumps)->a_loc );
    Free( (*mumps)->eltptr );
    Free( (*mumps)->eltvar );
    Free( (*mumps)->a_elt );
    Free( (*mumps)->perm_in );
    Free( (*mumps)->colsca );
    Free( (*mumps)->rowsca  );
    Free( (*mumps)->pivnul_list );
    Free( (*mumps)->listvar_schur );
    Free( (*mumps)->sym_perm );
    Free( (*mumps)->uns_perm );
    Free( (*mumps)->irhs_ptr);
    Free( (*mumps)->irhs_sparse);
    Free( (*mumps)->rhs_sparse);
    // Free( (*mumps)->rhs);
    Free( (*mumps)->redrhs);
    Free(*mumps);
  }
}
