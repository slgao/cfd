#include "boundary_val.h"
void boundaryvalues(
		    int imax,
		    int jmax,
		    int il,int ir,int jb,int jt,
		    double **U,
		    double **V
		    ){
  int i,j;
  /*bottom boundary*/
  if(jb==1){
    for(i=il-1; i<=ir; i++){
      V[i][0] = 0;
      U[i][0] = -U[i][1];
    }
   
  }
  /*top boundary*/
  if(jt==jmax){
    for(i=il-1; i<=ir; i++){
      V[i][jmax] = 0;
      U[i][jmax+1] = 2.0 -U[i][jmax];
    }
  }
  /*left boundary*/
  if(il==1){
    for(j=jb-1; j<=jt; j++){
      U[0][j] = 0;
      V[0][j] = -V[1][j];
    }
  }
  /*right boundary*/
  if(ir==imax){
    for(j=jb-1; j<=jt; j++){
      U[imax][j] = 0;
      V[imax+1][j] = -V[imax][j];
    }
  }
}





