#include "uvp.h"
#include "helper.h"
#include "math.h"

/* Computation of F and G */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
		  )
{
  double **duudx, **duvdy, **ddux, **dduy, **dvvdy, **duvdx, **ddvx, **ddvy;
  int i, j;
  duudx= matrix(0,imax-1,0,jmax); /* reserve memory for duudx .....*/
  duvdy= matrix(0,imax-1,0,jmax); /* change all matrix 1 to 0s */
  ddux= matrix(0,imax-1,0,jmax);
  dduy= matrix(0,imax-1,0,jmax);
  dvvdy= matrix(0,imax,0,jmax-1);
  duvdx= matrix(0,imax,0,jmax-1);
  ddvx= matrix(0,imax,0,jmax-1);
  ddvy= matrix(0,imax,0,jmax-1);

  for (i=1;i<imax;i++)
    for (j=1;j<jmax+1;j++)
      {
	duudx[i][j]=1/dx * ( pow( ( U[i][j]+U[i+1][j] )/2.0 , 2 ) - pow( ( U[i-1][j]+U[i][j] )/2.0 , 2 ) ) + 
	  alpha/dx * (fabs(U[i][j]+U[i+1][j]) * (U[i][j]-U[i+1][j])/4.0 - fabs(U[i-1][j]+U[i][j]) * (U[i-1][j]-U[i][j])/4.0 );
	duvdy[i][j]=1/dy * ( ( V[i][j]+V[i+1][j] ) * (U[i][j]+U[i][j+1] )/4.0 - ( V[i][j-1]+V[i+1][j-1] ) * ( U[i][j-1]+U[i][j] )/4.0 ) + 
	  alpha/dy * ( fabs( V[i][j]+V[i+1][j] ) * ( U[i][j]-U[i][j+1] )/4.0 - fabs( V[i][j-1]+V[i+1][j-1] ) * ( U[i][j-1]-U[i][j] )/4.0 );
	ddux[i][j]=( U[i+1][j]-2.0*U[i][j]+U[i-1][j] )/(dx * dx );
	dduy[i][j]=( U[i][j+1]-2.0*U[i][j]+U[i][j-1] )/(dy * dy );
      }
  for (i=1;i<imax+1;i++)
    for (j=1;j<jmax;j++)
      {
	dvvdy[i][j]=1/dy * ( pow( ( V[i][j]+V[i][j+1] )/2.0 , 2 ) - pow( ( V[i][j-1]+V[i][j] )/2.0 , 2 ) ) + 
	  alpha/dy * (fabs(V[i][j]+V[i][j+1]) * (V[i][j]-V[i][j+1])/4.0 - fabs(V[i][j-1]+V[i][j]) * (V[i][j-1]-V[i][j])/4.0 );
	duvdx[i][j]=1/dx * ( ( U[i][j]+U[i][j+1] ) * (V[i][j]+V[i+1][j] )/4.0 - ( U[i-1][j]+U[i-1][j+1] ) * ( V[i-1][j]+V[i][j] )/4.0 ) + 
	  alpha/dx * ( fabs( U[i][j]+U[i][j+1] ) * ( V[i][j]-V[i+1][j] )/4.0 - fabs( U[i-1][j]+U[i-1][j+1] ) * ( V[i-1][j]-V[i][j] )/4.0 );
	ddvx[i][j]=( V[i+1][j]-2.0*V[i][j]+V[i-1][j] )/(dx * dx );
	ddvy[i][j]=( V[i][j+1]-2.0*V[i][j]+V[i][j-1] )/(dy * dy );
      }

  for (i=1;i<imax;i++)   /* calculate F including the boundray value */
    for (j=1;j<jmax+1;j++)
      {
	F[i][j]=U[i][j] + dt * (1.0/Re * (ddux[i][j]+dduy[i][j]) - duudx[i][j] - duvdy[i][j] + GX);
      }
  for(j=1;j<jmax+1;j++)
    {
      F[0][j]=U[0][j];
      F[imax][j]=U[imax][j];
}     
                      
  for (i=1;i<imax+1;i++)   /* calculate G including the boundary value */
    for (j=1;j<jmax;j++)
      {
	G[i][j]=V[i][j] + dt * (1.0/Re * (ddvx[i][j]+ddvy[i][j]) - dvvdy[i][j] - duvdx[i][j] + GY);
      }             
  for(i=1;i<imax+1;i++)
    {
      G[i][0]=V[i][0];
      G[i][jmax]=V[i][jmax];
} 
}

/* calculate the right hand side of the pressure equation */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
		  )
{
  int i,j;
  for (i=1;i<imax+1;i++)   
    for (j=1;j<jmax+1;j++)
      {
	RS[i][j]=1.0/dt * ( ( F[i][j]-F[i-1][j] )/dx + (G[i][j]-G[i][j-1])/dy );
      }
}

/* calculate dt */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
		  )
{
  double a, b, c, umax, vmax, min;
  int i, j;
  a=Re/2.0 * 1.0/(1.0/(dx*dx) + 1.0/(dy*dy));
  umax=U[0][1];
  for (i=0;i<imax+1;i++) /* get the umax through the u values in the domain */ 
    for(j=1;j<jmax+1;j++)
      {
	if(fabs(U[i][j])>umax)
	  umax=fabs(U[i][j]);
      }
  b=dx/fabs(umax);

  vmax=V[0][1];
  for (i=1;i<imax+1;i++) /* get the vmax value through the v values in the domain */
    for(j=0;j<jmax+1;j++)
      {
	if(fabs(V[i][j])>vmax)
	  vmax=fabs(V[i][j]);
      }
  c=dy/fabs(vmax);

  if (a>=b)
    min=b;
  else
    min=a;
  if(min>=c)
    min=c;
  *dt=tau * min;
}

/* calculate U and V velocities next time step */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
		  )
{
  int i, j;
  for (i=1;i<imax;i++)
    for(j=1;j<jmax+1;j++)
      {
	U[i][j]=F[i][j]-dt/dx*(P[i+1][j]-P[i][j]);
	
} 
  for (i=1;i<imax+1;i++)
    for(j=1;j<jmax;j++)
      {
	V[i][j]=G[i][j]-dt/dy*(P[i][j+1]-P[i][j]);
      }
}
