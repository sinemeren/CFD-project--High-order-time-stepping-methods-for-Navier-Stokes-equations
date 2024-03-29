#ifndef __INIT_H_
#define __INIT_H_

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.hpp and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the 
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 */

#include <fstream>
int read_parameters( 
  std::string szFileName,
  double *Re,
  double *UI,
  double *VI,
  double *PI,
  double *GX,
  double *GY,
  double *t_end,
  double *xlength,
  double *ylength,
  double *dt,
  double *dx,
  double *dy,
  int  *imax,
  int  *jmax,
  double *alpha,
  double *omg,
  double *tau,
  int  *itermax,
  double *eps,
  double *dt_value);

void initCellInfo(const std::string &pgmFile, const int &imax, const int &jmax, matrix<unsigned int> &flag);
bool IsFluid(unsigned int val);
bool check_EW(int i, int j, matrix<int> geo_pgm);
bool chech_NS(int i, int j, matrix<int> geo_pgm);
void check_forbidcell(int imax, int jmax, matrix<int> geo_pgm);



#endif

