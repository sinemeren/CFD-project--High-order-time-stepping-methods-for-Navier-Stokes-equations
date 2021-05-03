#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax, int jmax, matrix<double> &U, matrix<double> &V );

void boundaryvalues_arbitrary(int imax, int jmax, matrix<double> &U, matrix<double> &V, const matrix<unsigned int> &flag);

int B_E(unsigned int flag);

int B_W(unsigned int flag);

int B_N(unsigned int flag);

int B_S(unsigned int flag);

int B_NE(unsigned int flag);

int B_NW(unsigned int flag);

int B_SE(unsigned int flag);

int B_SW(unsigned int flag);

#endif
