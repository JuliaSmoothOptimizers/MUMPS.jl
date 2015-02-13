#include "mpi.h"
#include "dmumps_c.h"
#include "zmumps_c.h"
#include "complex.h"

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

void* mumps_initialize_complex(int, int*, double*);
void  mumps_associate_matrix_complex(void*, int, int, double complex*, int*, int*);
void  mumps_factorize_complex(void*);
void  mumps_associate_rhs_complex(void*, int, double complex*);
void  mumps_solve_complex(void*, int*);
int   mumps_get_nrhs_complex(void*);
void  mumps_get_solution_complex(void*, double complex*);
void  mumps_get_info_complex(void*, int*, double*);
void  mumps_finalize_complex(void*);
void  mumps_alloc_complex(ZMUMPS_STRUC_C**);
void  mumps_free_complex(ZMUMPS_STRUC_C**);
