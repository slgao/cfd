#include "sor.h"
#include <math.h>
#include "Definitions.h"
#include "string.h"
#include "stdio.h"

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag,
  char *problem,
  double dp
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  int numcell;
  numcell=0;
  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      if(Flag[i][j]>=81 && Flag[i][j]<162)
	{
	  P[i][j] = (1.0-omg)*P[i][j]
	    + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
	  numcell++;
	}
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      if(Flag[i][j]>=81 && Flag[i][j]<162)
	{
	  rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
	    ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
	}
    }
  }
  rloc = rloc/numcell;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  /* set top and bottom boundary values */
  for(i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];
    P[i][jmax+1] = P[i][jmax];
  }
  if(strcmp(problem,"problem1.dat")==0)
    {
      for(j = 1; j <= jmax; j++) 
	{
	  P[0][j] = P[1][j];
	  P[imax+1][j] = P[imax][j];
	}
    }
  if(strcmp(problem,"problem3.dat")==0)
    {
      for(j = 1; j <= jmax; j++) 
	{
	  P[0][j] = P[1][j];
	  P[imax+1][j] = P[imax][j];
	}
    }
 if(strcmp(problem,"problem2.dat")==0)
   {
     for(j = 1; j <= jmax; j++)
       {
	 P[0][j] = 2*dp - P[1][j];
	 P[imax+1][j] = -P[imax][j];

	 /* P[0][j] = -P[1][j]; */
	 /* P[imax+1][j] =2*dp- P[imax][j]; */
       }
   }
 if(strcmp(problem,"test.dat")==0)
   {
     for(j = 1; j <= jmax; j++)
       {
 	 P[0][j] = P[1][j];
 	 P[imax+1][j] = P[imax][j];
       }
   }
 /* inner boundary */
  for (i=1;i<imax+1;i++)
    for (j=1;j<jmax+1;j++)
      {
	if(Flag[i][j]<81)
	  {
	    switch(Flag[i][j])
	      {
	      case B_N:   
		P[i][j]=P[i][j+1];
		break;
	      case B_O:   
		P[i][j]=P[i+1][j];
		break;
	      case B_S:   
		P[i][j]=P[i][j-1];
		break;
	      case B_W:   
		P[i][j]=P[i-1][j];
		break;
	      case B_NO:  
		P[i][j]=(P[i][j+1]+P[i+1][j])/2;
		break;
	      case B_SO:  
		P[i][j]=(P[i+1][j]+P[i][j-1])/2;
		break;
	      case B_SW:  
		P[i][j]=(P[i-1][j]+P[i][j-1])/2;
		break;
	      case B_NW:  
		P[i][j]=(P[i][j+1]+P[i-1][j])/2;
		break;
	      default:
	      /* 	printf("forbidden boundary cells exit!\n"); */
	      	break;
	      }
	  }
	else if(Flag[i][j]>=162)
	  {
	    switch(Flag[i][j])
	      {
	      case H_N:   
		P[i][j]=P[i][j+1];
		break;
	      case H_O:   
		P[i][j]=P[i+1][j];
		break;
	      case H_S:   
		P[i][j]=P[i][j-1];
		break;
	      case H_W:   
		P[i][j]=P[i-1][j];
		break;
	      case H_NO:  
		P[i][j]=(P[i][j+1]+P[i+1][j])/2;
		break;
	      case H_SO:  
		P[i][j]=(P[i+1][j]+P[i][j-1])/2;
		break;
	      case H_SW:  
		P[i][j]=(P[i-1][j]+P[i][j-1])/2;
		break;
	      case H_NW:  
		P[i][j]=(P[i][j+1]+P[i-1][j])/2;
		break;
	      default:
	      /* 	printf("forbidden boundary cells exit!\n"); */
	      	break;
	      }
	  }
      }
}

