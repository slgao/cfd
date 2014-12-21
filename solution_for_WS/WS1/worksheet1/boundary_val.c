#include "boundary_val.h"

/* set boundary values for U and V */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
		    )
{
  int i,j;
  for (i=1;i<imax+1;i++){
    V[i][0]=0;
    V[i][jmax]=0;
  }
  for (j=1;j<jmax+1;j++){
    U[0][j]=0;
    U[imax][j]=0;
  }

  for (i=1;i<imax+1;i++){
    U[i][0]=-U[i][1];
    U[i][jmax+1]=2-U[i][jmax];
}
  for (j=1;j<jmax+1;j++){
    V[0][j]=-V[1][j];
    V[imax+1][j]=-V[imax][j];
  }
}
