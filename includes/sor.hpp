#ifndef __SOR_H_
#define __SOR_H_

#include "datastructures.hpp"

/**
 * One GS iteration for the pressure Poisson equation. Besides, the routine must 
 * also set the boundary values for P according to the specification. The 
 * residual for the termination criteria has to be stored in res.
 * 
 * An \omega = 1 GS - implementation is given within sor.c.
 */
void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        matrix<double> &P,
        matrix<double> &RS,
        double *res
);

void sor_arbitrary(double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        matrix<double> &P,
        matrix<double> &RS,
        double *res,
		const matrix<unsigned int> flag);


#endif
