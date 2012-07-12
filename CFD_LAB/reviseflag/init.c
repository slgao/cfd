#include "helper.h"
#include "init.h"
#include "Definitions.h"

int read_parameters( char *argv[],       /* name of the file */
		     double *Re,                /* reynolds number   */
		     double *UI,                /* velocity x-direction */
		     double *VI,                /* velocity y-direction */
		     double *PI,                /* pressure */
		     double *GX,                /* gravitation x-direction */
		     double *GY,                /* gravitation y-direction */
		     double *t_end,             /* end time */
		     double *xlength,           /* length of the domain x-dir.*/
		     double *ylength,           /* length of the domain y-dir.*/
		     double *dt,                /* time step */
		     double *dx,                /* length of a cell x-dir. */
		     double *dy,                /* length of a cell y-dir. */
		     double *TI,                /* initial value of temperature */
		     double *Pr,                /* Prandtl number */
		     double *beta,              /* coefficient of thermal expansion */
		     int  *imax,                /* number of cells x-direction*/
		     int  *jmax,                /* number of cells y-direction*/
		     double *alpha,             /* uppwind differencing factor*/
		     double *omg,               /* relaxation factor */
		     double *tau,               /* safety factor for time step*/
		     int  *itermax,             /* max. number of iterations  */
		                                /* for pressure per time step */
		     double *eps,               /* accuracy bound for pressure*/
		     double *dt_value,          /* time for output */
		     
		     int *wl,
		     int *wr,
		     int *wt,
		     int *wb,
		     double *dp,
		     double *T_H,
		     double *T_C,
		     double *T_I)          
{
   READ_DOUBLE( argv[1], *xlength );
   READ_DOUBLE( argv[1], *ylength );

   READ_DOUBLE( argv[1], *Re    );
   READ_DOUBLE( argv[1], *t_end );
   READ_DOUBLE( argv[1], *dt    );

   READ_DOUBLE( argv[1], *TI    );/* for the energy properties */
   READ_DOUBLE( argv[1], *Pr    );
   READ_DOUBLE( argv[1], *beta    );
   READ_DOUBLE( argv[1], *T_H    );
   READ_DOUBLE( argv[1], *T_C    );
   READ_DOUBLE( argv[1], *T_I    );

   READ_INT   ( argv[1], *imax );
   READ_INT   ( argv[1], *jmax );

   READ_DOUBLE( argv[1], *omg   );
   READ_DOUBLE( argv[1], *eps   );
   READ_DOUBLE( argv[1], *tau   );
   READ_DOUBLE( argv[1], *alpha );

   READ_INT   ( argv[1], *itermax );
   READ_DOUBLE( argv[1], *dt_value );

   READ_DOUBLE( argv[1], *UI );
   READ_DOUBLE( argv[1], *VI );
   READ_DOUBLE( argv[1], *GX );
   READ_DOUBLE( argv[1], *GY );
   READ_DOUBLE( argv[1], *PI );
   /*READ_STRING( argv[1], problem);*/
   READ_INT   ( argv[1], *wl );
   READ_INT   ( argv[1], *wr );
   READ_INT   ( argv[1], *wt );
   READ_INT   ( argv[1], *wb );
   READ_DOUBLE( argv[1], *dp );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}

/* initialization of U V P */

void init_uvp(
  double UI,
  double VI,
  double PI,
  double TI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  double **TEMP,
  char *problem,
  int **Flag
	      )
{ 

  /* int i,j; */

  if(strcmp(problem,"problem3.dat")==0)
    {
      init_matrix(U,0,imax,1,jmax/2,0.0);
      init_matrix(U,0,imax,jmax/2+1,jmax,UI);
      init_matrix(V,1,imax,0,jmax,VI);
      init_matrix(P,1,imax,1,jmax,PI);
      init_matrix(TEMP,0,imax+1,0,jmax+1,TI);
    }
  else 
    {
      init_matrix(U,0,imax,1,jmax,UI);
      init_matrix(V,1,imax,0,jmax,VI);
      init_matrix(P,1,imax,1,jmax,PI);
      init_matrix(TEMP,0,imax+1,0,jmax+1,TI);
    }
}

void init_flag(char* problem, int imax, int jmax, int **Flag)
{
  int i,j;
  int **pic;
  pic=read_pgm(problem);
  for(i=0; i<=imax+1; i++)
    for(j=0; j<=jmax+1; j++)
      {
	/* if (pic[i][j]!=0)	 */
	/*   pic[i][j]=1; */
	/* else */
	/*   pic[i][j]=C_B; */
	if (pic[i][j]==0)
	  pic[i][j]=C_H;
	else if (pic[i][j]!=0 && pic[i][j]!=255)
	  pic[i][j]=C_B;
	else if (pic[i][j]==255)
	  pic[i][j]=1;
      }
    for (i=1; i<=imax; i++)
      for (j=1;j<=jmax;j++)
	{
	  /* Flag[i][j]=16*pic[i][j]+8*pic[i+1][j]+4*pic[i-1][j]+2*pic[i][j-1]+pic[i][j+1]; */
	  Flag[i][j]=81*pic[i][j]+27*pic[i+1][j]+9*pic[i-1][j]+3*pic[i][j-1]+pic[i][j+1];
	}
  
  
}






