#include "uvp.h"
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void calculate_dt(
		  double Re,
		  double tau,
		  double *dt,
		  double dx,
		  double dy,
		  int il,int ir,int jb,int jt,
		  double **U,
		  double **V
		  ){
  int i,j,myrank;
  double temp1;
  double temp2;
  double temp3;
  double u_tmp;
  double v_tmp;
  double u_max = fabs(U[il][jb]);
  double v_max = fabs(V[il][jb]);

  for(i=il-1; i <= ir; i++){
    for(j=jb; j <=jt; j++){
      u_max = fmax(fabs(U[i][j]), u_max);
    }
  }

  for(i=il; i <= ir; i++){
    for(j=jb-1; j <=jt; j++){
      v_max = fmax(fabs(V[i][j]), v_max);
    }
  }

  MPI_Reduce(&u_max,&u_tmp,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&v_max,&v_tmp,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0){
    u_max = u_tmp;
    v_max = v_tmp;
    temp1 = Re/2.0 * 1.0/(1.0/(dx*dx) + 1.0/(dy*dy));
    temp2 = dx/u_max;
    temp3 = dy/v_max;
    /* printf("temp1:%f,  temp2:%f,  temp3:%f\n",temp1,temp2,temp3); */
    *dt = tau * fmin(temp1, fmin(temp2, temp3));
  }
  MPI_Bcast(dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void calculate_fg(
		  double Re,
		  double GX,
		  double GY,
		  double alpha,
		  double dt,
		  double dx,
		  double dy,
		  int il,int ir,int jb,int jt,int imax,int jmax,
		  double **U,
		  double **V,
		  double **F,
		  double **G
		  ){
  int i,j;
  double u_2nd_deriv;
  double u_square_deriv;
  double u_uv_deriv;
  for(i=il-1; i<= ir; i++){
    for(j=jb; j<= jt; j++){
      u_2nd_deriv = 1.0/(dx*dx) * ( U[i+1][j] + U[i-1][j] - 2.0*U[i][j] )
	+ 1.0/(dy*dy) * ( U[i][j+1] + U[i][j-1] -2.0*U[i][j] );
      u_square_deriv = 1.0/dx * ( pow(0.5*(U[i][j]+U[i+1][j]),2.0)
				  - pow(0.5*(U[i-1][j]+U[i][j]),2.0) )
	+ alpha*1.0/dx * ( 0.5 * fabs(U[i][j]+U[i+1][j]) * 0.5 * (U[i][j]-U[i+1][j])
			   - 0.5 * fabs(U[i-1][j]+U[i][j]) * 0.5 * (U[i-1][j] - U[i][j]) );
      u_uv_deriv = 1.0/dy * ( 0.5 * (V[i][j]+V[i+1][j]) * 0.5 * (U[i][j]+U[i][j+1])
			      - 0.5 * (V[i][j-1]+V[i+1][j-1]) * 0.5 * (U[i][j-1]+U[i][j]) )
	+ alpha*1.0/dy * ( 0.5 * fabs(V[i][j]+V[i+1][j]) * 0.5 * (U[i][j]-U[i][j+1])
			   - 0.5 * fabs(V[i][j-1]+V[i+1][j-1]) * 0.5 * (U[i][j-1]-U[i][j]) );
      F[i][j] = U[i][j] + dt*( 1.0/Re*u_2nd_deriv - u_square_deriv - u_uv_deriv + GX );
    }
  }

  for(i=il; i<= ir; i++){
    for(j=jb-1; j<= jt; j++){
      double v_2nd_deriv = 1.0/(dx*dx) * ( V[i+1][j] + V[i-1][j] - 2.0*V[i][j] )
	+ 1.0/(dy*dy) * ( V[i][j+1] + V[i][j-1] -2.0*V[i][j] );
      double v_square_deriv = 1.0/dy * ( pow(0.5*(V[i][j]+V[i][j+1]),2.0)
					 - pow(0.5*(V[i][j-1]+V[i][j]),2.0) )
	+ alpha*1.0/dy * ( 0.5 * fabs(V[i][j]+V[i][j+1]) * 0.5 * (V[i][j]-V[i][j+1])
			   - 0.5 * fabs(V[i][j-1]+V[i][j]) * 0.5 * (V[i][j-1] - V[i][j]) );
      double v_uv_deriv = 1.0/dx * ( 0.5 * (V[i][j]+V[i+1][j]) * 0.5 * (U[i][j]+U[i][j+1])
				     - 0.5 * (V[i-1][j]+V[i][j]) * 0.5 * (U[i-1][j]+U[i-1][j+1]) )
	+ alpha*1.0/dx * ( 0.5 * fabs(U[i][j]+U[i][j+1]) * 0.5 * (V[i][j]-V[i+1][j])
			   - 0.5 * fabs(U[i-1][j]+U[i-1][j+1]) * 0.5 * (V[i-1][j]-V[i][j]) );
      G[i][j] = V[i][j] + dt*( 1.0/Re*v_2nd_deriv - v_square_deriv - v_uv_deriv + GY );
    }
  }
  /*set boundary value of F,G*/
  if(jb==1){
    for(i=il; i<= ir; i++){
      G[i][0] = V[i][0];
    }
  }
  if(jt==jmax){
    for(i=il; i<= ir; i++){
      G[i][jmax] = V[i][jmax];
    }
  }
  if(il==1){
    for(j=jb; j<= jt; j++){
      F[0][j] = U[0][j];
    }
  }
  if(ir==imax){
    for(j=jb; j<= jt; j++){
      F[imax][j] = U[imax][j];
    }
  }
}

void calculate_rs(
		  double dt,
		  double dx,
		  double dy,
		  int il,
		  int ir,int jb,int jt,
		  double **F,
		  double **G,
		  double **RS
		  ){
  int i,j;
  for(i=il; i<= ir; i++){
    for(j=jb; j<= jt; j++){
      RS[i][j] = 1.0/dt * ( 1.0/dx*(F[i][j]-F[i-1][j]) + 1.0/dy*(G[i][j]-G[i][j-1]) );
    }
  }
}

void calculate_uv(
		  double dt,
		  double dx,
		  double dy,
		  int il,int ir,int jb,int jt,
		  double **U,
		  double **V,
		  double **F,
		  double **G,
		  double **P
		  ){
  int i,j;
  for(i=il-1; i<= ir; i++){
    for(j=jb; j<= jt; j++){
      U[i][j] = F[i][j] - dt/dx* ( P[i+1][j] -P[i][j] );
    }
  }

  for(i=il; i<= ir; i++){
    for(j=jb-1; j<= jt; j++){
      V[i][j] = G[i][j] - dt/dy* ( P[i][j+1] - P[i][j] );
    }
  }
}





