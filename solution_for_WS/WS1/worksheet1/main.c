#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <math.h>
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
  double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, res, t,  **U, **V, **P, **RS, **F, **G;
  int imax, jmax, itermax, n, it;
  const char* dataFile = "cavity100.dat", * szProblem = "Visualize" ;
  
  /*Reading the data from .dat file*/
  read_parameters(dataFile,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value);
  
  /*Memory Allocation*/
  U = matrix(0,imax,0,jmax+1); 
  V = matrix(0,imax+1,0,jmax); 
  P = matrix(0,imax+1,0,jmax+1);
  RS = matrix(0,imax,0,jmax);
  F = matrix(0,imax,0,jmax); 
  G = matrix(0,imax,0,jmax);

  /* set t=0, n=0*/
  t=0.0;n=0;                             
  
  
  /* Initialization */
  init_uvp(UI,VI,PI,imax,jmax,U,V,P);
  
  /* Main Time loop */
  while(t<t_end)
    {
      /* Calculating dt */
      if (tau>0)
	calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
                  
      /* Boundary Conditions */
      boundaryvalues(imax,jmax,U,V);
      
      /* Calculating F and G */
      calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);

      /* Calculating Right Hand side */
      calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);


      /* set it and res */
      it = 0;                
      res = 15;

          	  
      /* loop for SOR iterations */
      while(it < itermax && res > eps)
	{
	  sor(omg,dx,dy,imax,jmax,P,RS,&res);
	  it++;
	}
      /* Calculating U and V from abouve calculated P*/
      calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
      /* output for visuialization */
      write_vtkFile(szProblem,n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
      /* Updating time and no of iterations*/
      t=t+dt;
      n=n+1;
    }
  

  /* free memory */
  free_matrix(U,0,imax,0,jmax+1);  
  free_matrix(V,0,imax+1,0,jmax); 
  free_matrix(P,0,imax+1,0,jmax+1);
  free_matrix(RS,0,imax,0,jmax);
  free_matrix(F,0,imax,0,jmax);
  free_matrix(G,0,imax,0,jmax);

  return -1;
}
