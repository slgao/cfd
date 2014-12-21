#ifndef __BOUNDARY_VAL_H__
#define __BOUNDARY_VAL_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int il,int ir,int jb,int jt,
  double **U,
  double **V
);

#endif

