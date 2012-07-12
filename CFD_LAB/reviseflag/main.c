#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <math.h>
#include <stdio.h>
#include "string.h"


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
int main(int argn, char* argv[]){
  double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, res, t, TI, Pr, beta, T_H, T_C,T_I,
    **U, **V, **P, **RS, **F, **G, **TEMP, **H, dp;
  int imax, jmax, itermax, n, it, wl, wr, wt ,wb, **Flag;
  /*  const char* dataFile = "cavity100.dat";*/
  const char * szProblem = "Visualize";
 
    char *p1 = "case11.pgm";

    char *p2= "case2.pgm";

    char *p3= "case3.pgm";

    char *p4 = "houseview.pgm";

  /*Reading the data from .dat file*/
    read_parameters(argv,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,&dx,&dy,&TI,&Pr,&beta,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value, &wl, &wr, &wt, &wb, &dp, &T_H, &T_C, &T_I);

  /*Memory Allocation*/
  Flag=imatrix(0,imax+1, 0, jmax+1);
  U = matrix(0,imax,0,jmax+1); 
  V = matrix(0,imax+1,0,jmax); 
  P = matrix(0,imax+1,0,jmax+1);
  TEMP = matrix(0,imax+1,0,jmax+1);
  RS = matrix(0,imax,0,jmax);
  F = matrix(0,imax,0,jmax); 
  G = matrix(0,imax,0,jmax);
  H = matrix(0,imax+1,0,jmax+1);

  /* set t=0, n=0*/
  t=0.0;n=0;                             
  printf("------------%s\n",argv[1]);
  
  /* Initialization */
  init_uvp(UI,VI,PI,TI,imax,jmax,U,V,P,TEMP,argv[1],Flag);
  if(strcmp(argv[1],"problem1.dat")==0)
    init_flag(p1,imax, jmax, Flag);  
  if(strcmp(argv[1],"problem2.dat")==0)
    init_flag(p2,imax, jmax, Flag);  
  if(strcmp(argv[1],"problem3.dat")==0)
    init_flag(p3,imax, jmax, Flag); 
  if(strcmp(argv[1],"test.dat")==0)
    init_flag(p4,imax, jmax, Flag);   


  /* init_flag(argv[1],imax, jmax, Flag);   */
  /* Main Time loop */
  while(t<t_end)
    {
      /* Calculating dt */
      if (tau>0)
	calculate_dt(Re,tau,&dt,dx,dy,Pr,imax,jmax,U,V);
                  
      /* Boundary Conditions */
      boundaryvalues(imax,jmax,U,V,TEMP,wl,wr,wt,wb,Flag,T_I);
      spec_boundary_value(argv[1],imax,jmax,T_H,T_C,U,V,TEMP);
      
      /* Calculating T(n+1) */
      calculate_TEMP(U,V,TEMP,Flag,imax,jmax,dt,dx,dy,alpha,Re,Pr,T_I);
      printf("temp is ::%f\n",TEMP[20][20]);

      /* Calculating F and G */
      calculate_fg(Re,GX,GY,alpha,dt,dx,dy,beta,imax,jmax,U,V,F,G,TEMP,Flag);

      /* Calculating Right Hand side */
      calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,Flag);


      /* set it and res */
      it = 0;                
      res = eps+1;

          	  
      /* loop for SOR iterations */
      while(it < itermax && res > eps)
	{
	  sor(omg,dx,dy,imax,jmax,P,RS,&res,Flag,argv[1],dp);
	  it++;
	}
      /* Calculating U and V from abouve calculated P*/
      calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,Flag);
      /* compute the HEAT */
      COMP_HEAT(U,V,TEMP,H,Flag,Re,Pr,imax,jmax,dx,dy);
      /* output for visuialization */
      if (n%10==0)
	write_vtkFile(szProblem,n,xlength,ylength,imax,jmax,dx,dy,U,V,TEMP,P);

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
  free_imatrix(Flag,0,imax+1,0,jmax+1);
  free_matrix(TEMP,0,imax+1,0,jmax+1);
  free_matrix(H,0,imax+1,0,jmax+1);

  return -1;
}
