#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_

/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeDensity(const double *const currentCell, double *density);

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double *const currentCell, const double * const density,double *velocity);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq);

#endif

