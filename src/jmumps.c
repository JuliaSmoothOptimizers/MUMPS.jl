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

// Initialize
// Initialize MUMPS data structure, store default parameters.
void* mumps_initialize_double(int sym, int* icntl, double* cntl) {

  DMUMPS_STRUC_C* pmumps;
  int i;

#ifdef JMUMPS_DEBUG
  int ierr, taskid, np;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  printf("MPI process %d out of %d\n", taskid, np);
#endif

  // Initialize MUMPS.
  mumps_alloc_double(&pmumps);
  if (pmumps == NULL) return NULL;

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

void* mumps_initialize_complex(int sym, int* icntl, double* cntl) {

  ZMUMPS_STRUC_C* pmumps;
  int i;

#ifdef JMUMPS_DEBUG
  int ierr, taskid, np;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  printf("MPI process %d out of %d\n", taskid, np);
#endif

  // Initialize MUMPS.
  mumps_alloc_complex(&pmumps);
  if (pmumps == NULL) return NULL;

  pmumps->sym = sym;
  pmumps->job = JOB_INIT;
  pmumps->par = 1;
  pmumps->comm_fortran = USE_COMM_WORLD;

  pmumps->n   = -1;
  pmumps->nz  = -1;

  zmumps_c(pmumps);

  // Fill in default parameter values.
  for (i = 0; i < 40; i++)
    (pmumps->icntl)[i] = icntl[i];

  for (i = 0; i < 5; i++)
    (pmumps->cntl)[i] = cntl[i];

  return (void*)pmumps;
}


// Associate matrix
// Associate pointers to matrix in Mumps data structure.
// This is a separate function so users can call it
// only on the host if necessary.
void mumps_associate_matrix_double(void* jmumps, int n, int nz,
                                   double* vals, int* irow, int* jcol) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;

#ifdef JMUMPS_DEBUG
    printf("Associating matrix with MUMPS struct initialized at %p\n", pmumps);
#endif

  pmumps->n  = n;
  pmumps->nz = nz;

  // irow/jcol are 1-based in Julia, as in Fortran.
  pmumps->irn = irow;
  pmumps->jcn = jcol;
  pmumps->a   = vals;

  return;
}

void mumps_associate_matrix_complex(void* jmumps, int n, int nz,
                                    double complex* vals, int* irow, int* jcol) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;

#ifdef JMUMPS_DEBUG
    printf("Associating matrix with MUMPS struct initialized at %p\n", pmumps);
#endif

  pmumps->n  = n;
  pmumps->nz = nz;

  // irow/jcol are 1-based in Julia, as in Fortran.
  pmumps->irn = irow;
  pmumps->jcn = jcol;
  pmumps->a   = (mumps_double_complex*)vals;

  return;
}


// Factorize
// Factorize input matrix.
void mumps_factorize_double(void* jmumps) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;

  pmumps->job = JOB_ANALYZE_FACTORIZE;
  dmumps_c(pmumps);

#ifdef JMUMPS_DEBUG
    printf("MUMPS factorized a matrix of size %d with %d nonzeros\n",
           pmumps->n, pmumps->nz);
#endif

    return;
}

void mumps_factorize_complex(void* jmumps) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;

  pmumps->job = JOB_ANALYZE_FACTORIZE;
  zmumps_c(pmumps);

#ifdef JMUMPS_DEBUG
    printf("MUMPS factorized a matrix of size %d with %d nonzeros\n",
           pmumps->n, pmumps->nz);
#endif

    return;
}


// Associate RHS
// Associate pointer to right-hand side in Mumps data structure.
// This is a separate function so users can call it
// only on the host if necessary.
void mumps_associate_rhs_double(void* jmumps, int nrhs, double* rhs) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;

  pmumps->nrhs = nrhs;
  pmumps->lrhs = pmumps->n;
  pmumps->rhs  = rhs;  // Will be overwritten with the solution.
}

void mumps_associate_rhs_complex(void* jmumps, int nrhs, double complex* rhs) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;

  pmumps->nrhs = nrhs;
  pmumps->lrhs = pmumps->n;
  pmumps->rhs  = (mumps_double_complex*)rhs;  // Will be overwritten with the solution.
}

// Solve
// Solve a linear system using the factorization computed previously.
void mumps_solve_double(void* jmumps, int* transposed) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int taskid, ierr;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  if (taskid == 0) {

#ifdef JMUMPS_DEBUG
    printf("Solving linear system at %p\n", pmumps);
    printf("  nrhs=%d\n", pmumps->nrhs);
    printf("     n=%d\n", pmumps->n);
    printf("    nz=%d\n", pmumps->nz);
    printf("   irn=%p\n", pmumps->irn);
    printf("   jcn=%p\n", pmumps->jcn);
    printf("     a=%p\n", pmumps->a);
    int k;
    for (k = 0; k < pmumps->n; k++) {
      printf("  rhs[%d]=%8.1e\n", k, pmumps->rhs[k]);
    }
    for (k = 0; k < pmumps->nz; k++)
      printf("  irow[%d]=%d  jcol[%d]=%d  a[%d]=%8.1e\n",
             k, (pmumps->irn)[k], k, (pmumps->jcn)[k], k, (pmumps->a)[k]);
#endif

    if (transposed != 0)
      pmumps->icntl[8] = 0;
  }

  pmumps->job  = JOB_SOLVE;
  dmumps_c(pmumps);
  return;
}

void mumps_solve_complex(void* jmumps, int* transposed) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;
  int taskid, ierr;

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  if (taskid == 0) {

#ifdef JMUMPS_DEBUG
    printf("Solving linear system at %p\n", pmumps);
    printf("  nrhs=%d\n", pmumps->nrhs);
    printf("     n=%d\n", pmumps->n);
    printf("    nz=%d\n", pmumps->nz);
    printf("   irn=%p\n", pmumps->irn);
    printf("   jcn=%p\n", pmumps->jcn);
    printf("     a=%p\n", pmumps->a);
    int k;
    for (k = 0; k < pmumps->n; k++) {
      printf("  rhs[%d]=%8.1e + (%8.1e)i\n", k, pmumps->rhs[k].r, pmumps->rhs[k].i);
    }
    for (k = 0; k < pmumps->nz; k++)
      printf("  irow[%d]=%d  jcol[%d]=%d  a[%d]=%8.1e + (%8.1e)i\n",
             k, (pmumps->irn)[k], k, (pmumps->jcn)[k], k, pmumps->a[k].r, pmumps->a[k].i);
#endif

    if (transposed != 0)
      pmumps->icntl[8] = 0;
  }

  pmumps->job  = JOB_SOLVE;
  zmumps_c(pmumps);
  return;
}


// Get solution
// Retrieve solution from Mumps data structure.
// This is a separate function so users can call it
// only on the host if necessary.
void mumps_get_solution_double(void* jmumps, double* x) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  int i;

  for (i=0; i < pmumps->n * pmumps->nrhs; i++)
    x[i] = pmumps->rhs[i];

  return;
}

void mumps_get_solution_complex(void* jmumps, double complex* x) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;
  ZMUMPS_COMPLEX  val;
  int i;

  for (i=0; i < pmumps->n * pmumps->nrhs; i++) {
    val = pmumps->rhs[i];
    x[i] = val.r + val.i * I;
  }

  return;
}


// Get nrhs
// Obtain the number of rhs associated with Mumps data structure.
int mumps_get_nrhs_double(void* jmumps) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;
  return pmumps->nrhs;
}

int mumps_get_nrhs_complex(void* jmumps) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;
  return pmumps->nrhs;
}

// Get info
// Obtain analysis / factorization / solve integer and real info vectors.
void mumps_get_info_double(void* jmumps, int* infog, double* rinfog) {

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

void mumps_get_info_complex(void* jmumps, int* infog, double* rinfog) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;
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
void mumps_finalize_double(void* jmumps) {

  DMUMPS_STRUC_C* pmumps = (DMUMPS_STRUC_C*)jmumps;

  if (pmumps == NULL) return;

#ifdef JMUMPS_DEBUG
  int taskid = 0, mpi_finalized = 0;
  MPI_Finalized(&mpi_finalized);

  if (!mpi_finalized)
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  if (taskid == 0)
    printf("Terminating MUMPS struct at %p\n", pmumps);
#endif

  pmumps->job = JOB_END;
  dmumps_c(pmumps);
  mumps_free_double(&pmumps);

#ifdef JMUMPS_DEBUG
  if (taskid == 0)
    printf("MUMPS instance terminated\n");
#endif

  return;
}

void mumps_finalize_complex(void* jmumps) {

  ZMUMPS_STRUC_C* pmumps = (ZMUMPS_STRUC_C*)jmumps;

  if (pmumps == NULL) return;

#ifdef JMUMPS_DEBUG
  int taskid = 0, mpi_finalized = 0;
  MPI_Finalized(&mpi_finalized);

  if (!mpi_finalized)
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  if (taskid == 0)
    printf("Terminating MUMPS struct at %p\n", pmumps);
#endif

  pmumps->job = JOB_END;
  zmumps_c(pmumps);
  mumps_free_complex(&pmumps);

#ifdef JMUMPS_DEBUG
  if (taskid == 0)
    printf("MUMPS instance terminated\n");
#endif

  return;
}


// Helper functions, strongly inspired by the MUMPS MATLAB interface.

void mumps_alloc_double(DMUMPS_STRUC_C** mumps){

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

void mumps_alloc_complex(ZMUMPS_STRUC_C** mumps){

  (*mumps) = malloc(sizeof(ZMUMPS_STRUC_C));
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

void mumps_free_double(DMUMPS_STRUC_C** mumps){
  if (*mumps != NULL){
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
    Free( (*mumps)->redrhs);
    Free(*mumps);
  }
}

void mumps_free_complex(ZMUMPS_STRUC_C** mumps){
  if (*mumps != NULL){
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
    Free( (*mumps)->redrhs);
    Free(*mumps);
  }
}
