#include "mpi.h"

/* spezielle Routinen fuer paralleles Nast */

void Programm_Message(char *txt);

void Programm_Sync(char *txt);

void Programm_Stop(char *txt);

void  init_parallel(int iproc, int jproc, int imax, int jmax, int myrank, int *il, int *ir, int *jb, int *jt, int *rank_le, int *rank_ri, int *rank_bo, int *rank_to, int *omg_i, int *omg_j, int num_proc);


void pressure_comm(double **P, int il, int ir, int jb, int jt, int rank_le, int rank_ri, int rank_bo, int rank_to, double *BufSend, double *BufRecv, MPI_Status *status, int chunk);

void uv_comm(double **U, double **V, int il, int ir, int jb, int jt, int rank_le, int rank_ri, int rank_bo, int rank_to, double *BufSend, double *BufRecv, MPI_Status *status, int chunk);

