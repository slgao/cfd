#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
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
  int **Flag,
  double T_I
		    );

void spec_boundary_value(
  char* problem,
  int imax,
  int jmax,
  double T_H,
  double T_C,
  double **U,
  double **V,
  double **TEMP

			 );

#endif
