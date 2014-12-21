#include "sor.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
#include <mpi.h>
void sor(
	 double omg,
	 double dx,
	 double dy,
	 int imax,int jmax,int il,int ir,int jb,int jt,int myrank,
	 double **P,
	 double **RS,
	 double *res
	 ){
  int i,j;
  double p_temp1;
  double p_temp2;
  double res_tmp;
  *res = 0;
  for(i=il; i<= ir; i++){
    for(j=jb; j<= jt; j++){
      P[i][j] = (1-omg)*P[i][j] + omg/2.0/( 1/(dx*dx) + 1/(dy*dy) )
	* ( 1/(dx*dx)*(P[i+1][j]+P[i-1][j]) +
	    1/(dy*dy)*(P[i][j+1]+P[i][j-1]) - RS[i][j] );
    }
  }

 
  MPI_Barrier(MPI_COMM_WORLD);


  /*set boundary values of P*/
  if(il==1){
    for(j=jb; j<= jt; j++){
      P[0][j] = P[1][j];
    }
  }
  if(ir==imax){
    for(j=jb; j<= jt; j++){
      P[imax+1][j] = P[imax][j];
    }
  }
 
  if(jb==1){
    for(i=il-1; i<= ir+1; i++){
      P[i][0] = P[i][1];
    
    }
  }

  if(jt==jmax){
    for(i=il; i<= ir; i++){
      P[i][jmax+1] = P[i][jmax];
     
    }
  }
 
  /* calculate residual */
  for(i=il; i<= ir; i++){
    for(j=jb; j<= jt; j++){
      p_temp1 = 1/(dx*dx) * (P[i+1][j]+P[i-1][j]-2.0*P[i][j]);
      p_temp2 = 1/(dy*dy) * (P[i][j+1]+P[i][j-1]-2.0*P[i][j]);
      *res += pow( (p_temp1 + p_temp2 -RS[i][j]), 2.0 ) / (imax*jmax);
   
    }
  }
 
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(res,&res_tmp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(myrank==0){
    *res = res_tmp;
    *res = sqrt(*res);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

}

