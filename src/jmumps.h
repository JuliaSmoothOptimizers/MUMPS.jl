#include "mpi.h"
#include "dmumps_c.h"

int   mumps_initialize_mpi(void);
int   mumps_finalize_mpi(void);
void* mumps_initialize(int, int*, double*);
void  mumps_factorize(void*, int, int, double*, int*, int*);
void  mumps_solve(void*, int, double*);
void  mumps_get_info(void*, int*, double*);
void  mumps_finalize(void*);
void  mumps_alloc(DMUMPS_STRUC_C**);
void  mumps_free(DMUMPS_STRUC_C**);

