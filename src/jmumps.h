#include "mpi.h"
#include "dmumps_c.h"

void* mumps_initialize(int, int*, double*);
void  mumps_associate_matrix(void*, int, int, double*, int*, int*);
void  mumps_factorize(void*);
void  mumps_associate_rhs(void*, int, double*);
void  mumps_solve(void*, int*);
int   mumps_get_nrhs(void*);
void  mumps_get_solution(void*, double*);
void  mumps_get_info(void*, int*, double*);
void  mumps_finalize(void*);
void  mumps_alloc(DMUMPS_STRUC_C**);
void  mumps_free(DMUMPS_STRUC_C**);

