#ifndef __UVP_HPP__
#define __UVP_HPP__


#include "datastructures.hpp"
#include <string>
/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */

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
        const matrix<double> &U,
        const matrix<double> &V,
        matrix<double> &F,
        matrix<double> &G
);

void calculate_fg_arbitrary(
        double Re,
        double GX,
        double GY,
        double alpha,
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U,
		matrix<double> &V,
        matrix<double> &F,
        matrix<double> &G,
        const matrix<unsigned int> &flags
);



/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U_star,
        matrix<double> &V_star,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS
);

void calculate_rs_arbitrary(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
		matrix<double> &U_star,
		matrix<double> &V_star,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS,
        const matrix<unsigned int> &flags
);



/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(double Re,
        double tau,
        double *dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        const matrix<double> &U,
        const matrix<double> &V);

void calculate_dt_arbitrary(double Re,
        double tau,
        double *dt,
        double dx,
        double dy,
        int imax,
        int jmax,
		const matrix<double> &U,
		const matrix<double> &V);


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */

void calculate_uv(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U,
        matrix<double> &V,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        std::string* method
);

void calculate_uv_arbitrary(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U,
        matrix<double> &V,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        std::string* method,
		const matrix<unsigned int> &flags

);


void calculate_uv_star(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U,
        matrix<double> &V,
        matrix<double> &U_star,
        matrix<double> &V_star,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        std::string &method,
        int k
);

void calculate_uv_star_arbitrary(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &U,
        matrix<double> &V,
        matrix<double> &U_star,
        matrix<double> &V_star,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        std::string &method,
        int k,
		const matrix<unsigned int> &flags
);

void calculate_uv_k(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &P,
        matrix<double> &F,
        matrix<double> &G,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        int k
);

void calculate_uv_k_arbitrary(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &P,
        matrix<double> &F,
        matrix<double> &G,
        matrix3D<double> &UU,
        matrix3D<double> &VV,
        int k,
		const matrix<unsigned int> &flags
);

#endif
