#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "sor.hpp"
#include <cstdio>
#include <iostream>
#include "uvp.hpp"
#include "boundary_val.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - DONE - read the program configuration file using read_parameters()
 * - DONE - set up the matrices (arrays) needed. Use the predefined matrix<typename> type and give initial values in the constructor.
 * - DONE - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the 
 * operations, a definition is defined already within uvp.h):
 *
 * - DONE - calculate_dt() Determine the maximal time step size.
 * - DONE - boundaryvalues() Set the boundary values for the next time step.
 * - DONE - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - DONE - calculate_rs()
 * - DONE - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - DONE - calculate_uv() Calculate the velocity at the next time step.
 */


int main(int argn, char* args[]){

	std::string problem;
	std::string method;
	std::string fileNameTime, fileNameResidual, fileNameIter;

	if (argn == 3){

	    problem = args[1];
		method = args[2];

		if ( ((problem != "LidDrivenCavity") && (problem != "FlowOverStep") && (problem != "PlaneShearFlow") &&  (problem != "KarmanVortexStreet"))  &&  ((method != "ExplicitEuler") && (method != "Heun") && (method != "RungeKutta")) ){
			printf(
					"Incorrect usage! Please use the following format:"
					" ./sim [problem_name] [method] \n"
					"Available problems:"
					"\tLidDrivenCavity, "
					"FlowOverStep, "
					"PlaneShearFlow, "
					"KarmanVortexStreet, "
					"\nAvailable methods:"
					"\tExplicitEuler, "
					"Heun, "
					"RungeKutta"
					"\nExample usage:"
					"\t\t./sim LidDrivenCavity ExplicitEuler \n"
					);
			return EXIT_FAILURE;
		}

		else if ((problem != "LidDrivenCavity") && (problem != "FlowOverStep") && (problem != "PlaneShearFlow") && (problem != "KarmanVortexStreet")){
			printf(
					"Incorrect problem name! Please use one of the following problem names:"
					"\tLidDrivenCavity, "
					"FlowOverStep, "
					"PlaneShearFlow, "
					"KarmanVortexStreet, "
					"\nPlease use the following format:"
					"\t\t\t\t\t./sim [problem_name] [method] \n"
					"Example usage:"
					"\t\t\t\t\t\t\t\t./sim LidDrivenCavity ExplicitEuler \n"
					);
			return EXIT_FAILURE;

		}
		else if((method != "ExplicitEuler") && (method != "Heun") && (method != "RungeKutta")){
			printf(
					"Incorrect method! Please use one of the following methods:"
					"\tExplicitEuler, "
					"Heun, "
					"RungeKutta"
					"\nPlease use the following format:"
					"\t\t\t\t./sim [problem_name] [method]\n"
					"Example usage:"
					"\t\t\t\t\t\t\t./sim LidDrivenCavity ExplicitEuler \n "
					);
			return EXIT_FAILURE;
		}

	}
   else {
		printf(
				"Incorrect usage! Please use the following format:"
				" ./sim [problem_name] [method] \n"
				"Available problems:"
				"\tLidDrivenCavity, "
				"FlowOverStep, "
				"PlaneShearFlow, "
				"KarmanVortexStreet, "
				"\nAvailable methods:"
				"\tExplicitEuler, "
				"Heun, "
				"RungeKutta"
				"\nExample usage:"
				"\t\t./sim LidDrivenCavity ExplicitEuler \n"
				);
		return EXIT_FAILURE;
  }


        // File names
	fileNameTime = problem + "_" + method + "_Time.txt";
	fileNameResidual = problem + "_" + method + "_Residual.txt";
	fileNameIter = problem + "_" +  method + "_Iter.txt";
	
	// Folder Creation
	DIR* dir = opendir("../Results");
	if (dir) {
	    /* Directory exists. */
	    closedir(dir);
	} else if (ENOENT == errno) {
	    /* Directory does not exist. */
		mkdir("../Results", 0777);
	} else {
		printf("Results folder cannot be created!");
	}

	std::string dir_problem = "../Results/" + problem;
	DIR* dir1 = opendir(dir_problem.c_str());
	if (dir1) {
		/* Directory exists. */
		closedir(dir1);
	} else if (ENOENT == errno) {
		/* Directory does not exist. */
		mkdir(dir_problem.c_str(), 0777);
	} else {
		printf("Problem folder cannot be created!");
	}

	std::string dir_method = dir_problem + "/" + method;
	DIR* dir2 = opendir(dir_method.c_str());
	if (dir2) {
		/* Directory exists. */
		closedir(dir2);
	} else if (ENOENT == errno) {
		/* Directory does not exist. */
		mkdir(dir_method.c_str(), 0777);
	} else {
		printf("Method folder cannot be created!");
	}

	// decleration of Variables
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;

    int imax, jmax, itermax, numofk;
    std::string szFileName = "../Source/dat/" + problem +".dat";
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, & GY, &t_end, &xlength, &ylength, &dt, & dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);



    // setting up the matrices and assigning initial values
    matrix<double> F(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<double> G(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<double> RS(imax+2, std::vector<double>(jmax+2, 0.0));

    matrix<double> P(imax+2, std::vector<double>(jmax+2, PI));
    matrix<double> U(imax+2, std::vector<double>(jmax+2, UI));
    matrix<double> V(imax+2, std::vector<double>(jmax+2, VI));

    matrix<double> U_star(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<double> V_star(imax+2, std::vector<double>(jmax+2, 0.0));

    matrix<unsigned int> Flags(imax+2, std::vector<unsigned int>(jmax+2, 0));

    // setting up the k values
    if(method == "ExplicitEuler"){ numofk = 1; }
    else if(method == "Heun"){ numofk = 2;}
    else if (method == "RungeKutta") { numofk = 4; }

    matrix3D<double> UU(numofk, std::vector<std::vector<double>>(imax+2, std::vector<double>(jmax+2, 0.0)));
    matrix3D<double> VV(numofk, std::vector<std::vector<double>>(imax+2, std::vector<double>(jmax+2, 0.0)));

    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, & GY, &t_end, &xlength, &ylength, &dt, & dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);
    // reading parameters

	if (problem != "LidDrivenCavity"){
	std::string pgm_file = "../Source/pgm/" + problem + ".pgm";
	initCellInfo(pgm_file,imax,jmax, Flags);
	}

    // arrays to save SOR iterations, time, residual
    std::vector<double> SOR, time, residual;
    SOR.reserve(itermax); time.reserve(itermax); residual.reserve(itermax);

    // initializing t and iteration variable for writing VTK file
    double t = 0.0;
    int iter_visual = 0;
    int iter;
    double res;

	// performing iterations
	while (t <= t_end) {

		// calculation of maximal step size
		if (problem == "LidDrivenCavity"){
		//calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		dt = 0.005;
		} else{
		calculate_dt_arbitrary(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		dt = 0.05;
		}

		// setting the boundary values
		if (problem == "LidDrivenCavity"){
		boundaryvalues(imax, jmax, U, V);
		} else {
		boundaryvalues_arbitrary(imax, jmax, U, V, Flags);
		}

		for (size_t k = 0; k < numofk; ++k) {

			if (problem == "LidDrivenCavity"){
			calculate_uv_star(dt, dx, dy, imax, jmax, U, V, U_star, V_star, UU, VV, method, k);
			} else {
			calculate_uv_star_arbitrary(dt, dx, dy, imax, jmax, U, V, U_star, V_star, UU, VV, method, k, Flags);
			}

			if (problem == "LidDrivenCavity"){
			boundaryvalues(imax, jmax, U_star, V_star);
			} else {
			boundaryvalues_arbitrary(imax,jmax, U_star, V_star, Flags);
			}

			// calculation of F and G matrices
			if (problem == "LidDrivenCavity"){
			calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U_star, V_star, F, G);
			} else{
			calculate_fg_arbitrary(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, U_star, V_star, F, G, Flags);
			}

			// calculation of right hand side
			if (problem == "LidDrivenCavity"){
			calculate_rs(dt, dx, dy, imax, jmax,U_star, V_star, F, G, RS);
			} else {
			calculate_rs_arbitrary(dt, dx, dy, imax, jmax,U_star, V_star, F, G, RS, Flags);
			}


			// performing pressure poisson iteration
			iter = 0;

			res = 5.0;

			while ((iter < itermax) && (eps < res)) {

				if (problem == "LidDrivenCavity"){
				sor(omg, dx, dy, imax, jmax, P, RS, &res);
				} else {
				sor_arbitrary(omg, dx, dy, imax, jmax, P, RS, &res, Flags);
				}

				iter = iter + 1;

			}

			residual.push_back(res); time.push_back(t); SOR.push_back(iter);

			if (problem == "LidDrivenCavity"){
			calculate_uv_k(dt, dx, dy, imax, jmax, P, F, G, UU, VV, k);
			} else {
			calculate_uv_k_arbitrary(dt, dx, dy, imax, jmax, P, F, G, UU, VV, k, Flags);
			}

		}

		// calculation of velocities
		if (problem == "LidDrivenCavity"){
		calculate_uv(dt, dx, dy, imax, jmax, U, V, UU, VV, &method);
		} else {
		calculate_uv_arbitrary(dt, dx, dy, imax, jmax, U, V, UU, VV, &method, Flags);
		}

		//update for time
		t = t + dt;


		// writing VTK file
		bool IsDoubleEqual = fabs(t - (iter_visual * dt_value)) < 0.0001;
		if ((t != 0.0) && ((t > (iter_visual * dt_value)) || (IsDoubleEqual))) {

			write_vtkFile((dir_method + "/" + problem + "." + method).c_str(), iter_visual + 1, xlength, ylength, imax, jmax, dx, dy, U, V, P);
			iter_visual = iter_visual + 1;

		}

		printf("t = %f ,dt = %f, Res = %f,iterations=%d \n",t,dt,res,iter-1);

	}



	// outputs for plotting
	saveToFile(residual, fileNameResidual, problem, method);
	saveToFile(time, fileNameTime, problem, method);
	saveToFile(SOR, fileNameIter, problem, method);

    return -1;
}
