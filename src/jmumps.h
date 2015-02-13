#include "mpi.h"
#include "dmumps_c.h"

void* mumps_initialize_double(int, int*, double*);
void  mumps_associate_matrix_double(void*, int, int, double*, int*, int*);
void  mumps_factorize_double(void*);
void  mumps_associate_rhs_double(void*, int, double*);
void  mumps_solve_double(void*, int*);
int   mumps_get_nrhs_double(void*);
void  mumps_get_solution_double(void*, double*);
void  mumps_get_info_double(void*, int*, double*);
void  mumps_finalize_double(void*);
void  mumps_alloc_double(DMUMPS_STRUC_C**);
void  mumps_free_double(DMUMPS_STRUC_C**);

