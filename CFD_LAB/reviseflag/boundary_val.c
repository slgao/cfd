#include "boundary_val.h"
#include "helper.h"
#include "Definitions.h"
#include "string.h"
#include "stdio.h"

/* set boundary values for U and V */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  double **TEMP,
  int wl,
  int wr,
  int wt,
  int wb,
  int ** Flag,
  double T_I
		    )
{
  int i,j;
  switch(wl)
    {
    case 1:
      for (j=1;j<jmax+1;j++)
	{
	  U[0][j]=0.0;
	  V[0][j]=-V[1][j];
	}
      break;
    case 2:
      for (j=1;j<jmax+1;j++)
	{
	  V[0][j]=V[1][j];
	  U[0][j]=0.0;
	}
      break;
    case 3:
      for (j=1;j<jmax+1;j++)
	{
	  U[0][j]=U[1][j];
	  V[0][j]=V[1][j];
	}
      break;
    }
  switch(wr)
    {
    case 1:
      for (j=1;j<jmax+1;j++)
	{
	  V[imax+1][j]=-V[imax][j];
	  U[imax][j]=0.0;
	}
      break;
    case 2:
      for (j=1;j<jmax+1;j++)
	{
	  U[imax][j]=0.0;
	  V[imax+1][j]=V[imax][j];
	}
      break;
    case 3:
      for (j=1;j<jmax+1;j++)
	{
	  U[imax][j]=U[imax-1][j];
	  V[imax+1][j]=V[imax][j];
	}
      break;
    }
  switch(wt)
    {
    case 1:
      for (i=1;i<imax+1;i++)
	{
	  U[i][jmax+1]=-U[i][jmax];
	  V[i][jmax]=0.0;
	}
      break;
    case 2:
      for (i=1;i<imax+1;i++)
	{
	  U[i][jmax+1]=U[i][jmax];
	  V[i][jmax]=0.0;
	}
      break;
    case 3:
      for (i=1;i<imax+1;i++)
	{
	  U[i][jmax+1]=U[i][jmax];
	  V[i][jmax]=V[i][jmax-1];
	}
      break;
    }
  switch(wb)
    {
    case 1:
      for (i=1;i<imax+1;i++)
	{
	  U[i][0]=-U[i][1];
	  V[i][0]=0.0;
	}
      break;
    case 2:
      for (i=1;i<imax+1;i++)
	{
	  U[i][0]=U[i][1];
	  V[i][0]=0.0;
	}
      break;
    case 3:
      for (i=1;i<imax+1;i++)
	{
	  U[i][0]=U[i][1];
	  V[i][0]=V[i][1];
	}
      break;
    }
  /* set all walls to be adiabatic */
  for(j=1;j<=jmax;j++)
    {
      TEMP[0][j] = TEMP[1][j];  /* dT/dn = 0 */ 
      TEMP[imax+1][j] = TEMP[imax][j];
    }
  for (i=1;i<=imax;i++)
    {
      TEMP[i][0] = TEMP[i][1];
      TEMP[i][jmax+1] = TEMP[i][jmax]; 
    }

  /* setting boundaries for the inner boundary cells */
  for (i=1;i<imax+1;i++)
    for(j=1;j<jmax+1;j++)
      {
	if(Flag[i][j]<81) /* filter the obstacle cells */
	  {
	    switch (Flag[i][j])
	      {
	      case B_N:
		V[i][j]   = 0.0;
		U[i][j]   = -U[i][j+1];
		U[i-1][j] = -U[i-1][j+1];
		TEMP[i][j] = TEMP[i][j+1];
		break;
	      case B_O:
		U[i][j]   = 0.0;
		V[i][j]   = -V[i+1][j];
		V[i][j-1] = -V[i+1][j-1];
		TEMP[i][j] = TEMP[i+1][j];
		break;
	      case B_S:
		V[i][j-1] = 0.0;
		U[i][j]   = -U[i][j-1];
		U[i-1][j] = -U[i-1][j-1];
		TEMP[i][j] = TEMP[i][j-1];
		break;
	      case B_W:
		U[i-1][j] = 0.0;
		V[i][j]   = -V[i-1][j];
		V[i][j-1] = -V[i-1][j-1];
		TEMP[i][j] = TEMP[i-1][j];
		break;
	      case B_NO:
		V[i][j]   = 0.0;
		U[i][j]   = 0.0;
		V[i][j-1] = -V[i+1][j-1];
		U[i-1][j] = -U[i-1][j+1];
		TEMP[i][j] = 0.5*(TEMP[i][j+1]+TEMP[i+1][j]);
		break;
	      case B_SO:
		V[i][j-1] = 0.0;
		U[i][j]   = 0.0;
		V[i][j]   = -V[i+1][j];
		U[i-1][j] = -U[i-1][j-1];
		TEMP[i][j] = 0.5*(TEMP[i][j-1]+TEMP[i+1][j]);
		break;
	      case B_SW:
		V[i][j-1] = 0.0;
		U[i-1][j] = 0.0;
		V[i][j]   = -V[i-1][j];
		U[i][j]   = -U[i][j-1];
		TEMP[i][j] = 0.5*(TEMP[i][j-1]+TEMP[i-1][j]);
		break;
	      case B_NW:
		V[i][j]   = 0.0;
		U[i-1][j] = 0.0;
		V[i][j-1] = -V[i-1][j-1];
		U[i][j]   = -U[i][j+1];
		TEMP[i][j] = 0.5*(TEMP[i][j+1]+TEMP[i-1][j]);
		break;
	      default:
	      /* 	printf("forbidden boundary cells exits!\n"); */
	      	break;
	      }
	  }

	else if (Flag[i][j]>=162) /* filter the heating inner cells */
	  {
	    switch (Flag[i][j])
	      {
	      case H_N:
		V[i][j]   = 0.0;
		U[i][j]   = -U[i][j+1];
		U[i-1][j] = -U[i-1][j+1];
		TEMP[i][j] =2*T_I-TEMP[i][j+1];
		break;
	      case H_O:
		U[i][j]   = 0.0;
		V[i][j]   = -V[i+1][j];
		V[i][j-1] = -V[i+1][j-1];
		TEMP[i][j] =2*T_I-TEMP[i+1][j];
		break;
	      case H_S:
		V[i][j-1] = 0.0;
		U[i][j]   = -U[i][j-1];
		U[i-1][j] = -U[i-1][j-1];
		TEMP[i][j] = 2*T_I-TEMP[i][j-1];
		break;
	      case H_W:
		U[i-1][j] = 0.0;
		V[i][j]   = -V[i-1][j];
		V[i][j-1] = -V[i-1][j-1];
		TEMP[i][j] = 2*T_I-TEMP[i-1][j];
		break;
	      case H_NO:
		V[i][j]   = 0.0;
		U[i][j]   = 0.0;
		V[i][j-1] = -V[i+1][j-1];
		U[i-1][j] = -U[i-1][j+1];
		TEMP[i][j] = 2*T_I-0.5*(TEMP[i][j+1]+TEMP[i+1][j]);
		break;
	      case H_SO:
		V[i][j-1] = 0.0;
		U[i][j]   = 0.0;
		V[i][j]   = -V[i+1][j];
		U[i-1][j] = -U[i-1][j-1];
		TEMP[i][j] = 2*T_I-0.5*(TEMP[i][j-1]+TEMP[i+1][j]);
		break;
	      case H_SW:
		V[i][j-1] = 0.0;
		U[i-1][j] = 0.0;
		V[i][j]   = -V[i-1][j];
		U[i][j]   = -U[i][j-1];
		TEMP[i][j] = 2*T_I-0.5*(TEMP[i][j-1]+TEMP[i-1][j]);
		break;
	      case H_NW:
		V[i][j]   = 0.0;
		U[i-1][j] = 0.0;
		V[i][j-1] = -V[i-1][j-1];
		U[i][j]   = -U[i][j+1];
		TEMP[i][j] = 2*T_I-0.5*(TEMP[i][j+1]+TEMP[i-1][j]);
		break;
	      default:
	      /* 	printf("forbidden boundary cells exits!\n"); */
	      	break;
	      }
	  }

      }

}
void spec_boundary_value(
 char* problem,
 int imax,
 int jmax,
 double T_H,
 double T_C,
 double **U,
 double **V,
 double **TEMP

			 )
{
  if (strcmp(problem, "problem1.dat")==0)
    {
      int j;
      for(j=1;j<jmax+1;j++)
        {
	  U[0][j]=1.0;
	}
      for(j=0;j<jmax+1;j++)
	{
	  V[0][j]=-V[1][j];
	}
    }
  else if (strcmp(problem, "problem2.dat")==0)
    {
      int j;
      for(j=1;j<jmax+1;j++)
      	U[0][j]=1.0;
      for(j=0;j<=jmax;j++)
      	V[0][j]=-V[1][j];
    }
  else if (strcmp(problem, "problem3.dat")==0)
    {
      int j;
      for(j=jmax/2+1;j<jmax+1;j++)
        {
	  U[0][j]=1.0;
	}
      for(j=jmax/2;j<=jmax;j++)
	{
	  V[0][j]=-V[1][j];
	}
    }
  else if(strcmp(problem, "test.dat")==0)
    {
      int j;
      int i;
     for(j=1; j<=jmax; j++)
       {
	TEMP[0][j] =  TEMP[1][j];         
	TEMP[imax+1][j] = TEMP[imax][j]; /* adiabatic walls */
       }
     for(i=1;i<=imax;i++)
       {
	TEMP[i][0] =2*T_H-TEMP[i][1];          /* lower wall heated */
	TEMP[i][jmax+1] =2*T_C - TEMP[i][jmax]; /* upper wall cooled */
       }
    
    }

}














