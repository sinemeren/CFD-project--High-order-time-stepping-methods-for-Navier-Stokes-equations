#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include "boundary_val.hpp"

// Determines the value of F and G
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
        matrix<double> &G)
{

    for( size_t i = 1; i <= imax-1; i++){

        for( size_t j = 1; j <= jmax; j++){

            F.at(i).at(j) =

            				(1.0 / Re) *( ( ( U.at(i-1).at(j) - 2.0 * U.at(i).at(j) + U.at(i+1).at(j) ) / ( dx * dx ) ) +  ( ( U.at(i).at(j-1) - 2.0 * U.at(i).at(j) + U.at(i).at(j+1) ) / ( dy * dy ) ) )

							// d(u^2)/d(x)
							- (0.25) * ( 1.0 / dx ) * (  pow( ( U.at(i).at(j) + U.at(i+1).at(j) ), 2.0 )  -   pow( ( U.at(i-1).at(j) + U.at(i).at(j) ), 2.0 )
                            + alpha * ( ( fabs( U.at(i).at(j) + U.at(i+1).at(j) )  * ( U.at(i).at(j) - U.at(i+1).at(j)) ) -   (  fabs( U.at(i-1).at(j) + U.at(i).at(j)) * ( U.at(i-1).at(j) - U.at(i).at(j) ) ) ) )

							// d(uv)/d(y)
							- (0.25) * ( 1.0 / dy ) * ( ( V.at(i).at(j) + V.at(i+1).at(j) ) * ( U.at(i).at(j) + U.at(i).at(j+1) ) - ( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) + U.at(i).at(j) )
                            + alpha * ( fabs(V.at(i).at(j) + V.at(i+1).at(j)) * ( U.at(i).at(j) - U.at(i).at(j+1) ) - fabs( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) - U.at(i).at(j) ) ) )
                            + GX;

        }

    }



    for(size_t i = 1; i <= imax; i++){

        for( size_t j = 1; j <= jmax - 1; j++){

            G.at(i).at(j) =

							(1.0 / Re) *( ( ( V.at(i-1).at(j) - 2.0 * V.at(i).at(j) + V.at(i+1).at(j) ) / ( dx * dx ) ) +  ( ( V.at(i).at(j-1) - 2.0 * V.at(i).at(j) + V.at(i).at(j+1) ) / ( dy * dy ) ) )

							// d(v^2) / dy
							- (0.25) * ( 1.0 / dy ) * (  pow( ( V.at(i).at(j) + V.at(i).at(j+1) ), 2.0 )  -   pow( ( V.at(i).at(j-1) + V.at(i).at(j) ), 2.0 )
														 + alpha * ( ( fabs( V.at(i).at(j) + V.at(i).at(j+1) )  * ( V.at(i).at(j) - V.at(i).at(j+1)) ) -   (  fabs( V.at(i).at(j-1) + V.at(i).at(j)) * ( V.at(i).at(j-1) - V.at(i).at(j) ) ) ) )

							// d(uv)/dx
							- (0.25) * ( 1.0 / dx ) * ( ( U.at(i).at(j) + U.at(i).at(j+1) ) * ( V.at(i).at(j) + V.at(i+1).at(j) ) - ( U.at(i-1).at(j) + U.at(i-1).at(j+1)) * ( V.at(i-1).at(j) + V.at(i).at(j) )
							+ alpha * ( fabs(U.at(i).at(j) + U.at(i).at(j+1)) * ( V.at(i).at(j) - V.at(i+1).at(j) ) - fabs( U.at(i-1).at(j) + U.at(i-1).at(j+1)) * ( V.at(i-1).at(j) - V.at(i).at(j) ) ) )
							+ GY;

        }

    }


    for (int i = 1; i <=imax ; i++) {
    	//Boundary conditions for G
        G.at(i).at(0) = V.at(i).at(0);
        G.at(i).at(jmax) = V.at(i).at(jmax);

        }

        //Boundary conditions for F
        for (int j = 1; j <=jmax ; j++) {
        F.at(0).at(j) = U.at(0).at(j);
        F.at(imax).at(j) = U.at(imax).at(j);
        }



}

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
){

	/* ---------------------------------------set boundary values for F and G-------------------------------------------------------*/

	// -------------------------------------------------- Left boundary --------------------------------------------------------------
	for(size_t j = 1; j <= jmax; j++){

		if( B_N(flags.at(0).at(j)) ){ G.at(0).at(j) = V.at(0).at(j); }

		if( B_S(flags.at(0).at(j)) ){ G.at(0).at(j-1) = V.at(0).at(j-1); }

		if( B_W(flags.at(0).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_W at left boundary in FG"));}

		if( B_E(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); }

		if( B_NE(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); G.at(0).at(j) = V.at(0).at(j); }

		if( B_NW(flags.at(0).at(j)) ){throw std::runtime_error(std::string("impossible to have B_NW at left boundary in FG")); }

		if( B_SE(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); G.at(0).at(j-1) = V.at(0).at(j-1); }

		if( B_SW(flags.at(0).at(j))){throw std::runtime_error(std::string("impossible to have B_SW at left boundary in FG")); }

		// Inflow - only x direction
		if((flags.at(0).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(j) = U.at(0).at(j);}

		// Outflow - only x direction
		if((flags.at(0).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at left boundary in FG")); }

	}

	// ---------------------------------------------------- Right boundary -------------------------------------------------------------
	for(size_t j = 1; j <= jmax; j++){

		if( B_N(flags.at(imax+1).at(j)) ){ G.at(imax+1).at(j) = V.at(imax+1).at(j); }

		if( B_S(flags.at(imax+1).at(j)) ){ G.at(imax+1).at(j-1) = V.at(imax+1).at(j-1); }

		if( B_W(flags.at(imax+1).at(j)) ){ F.at(imax).at(j) = U.at(imax).at(j); }

		if( B_E(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_E at right boundary in FG")); }

		if( B_NE(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_NE at right boundary in FG")); }

		if( B_NW(flags.at(imax+1).at(j)) ){ F.at(imax).at(j) = U.at(imax).at(j); G.at(imax+1).at(j) = V.at(imax+1).at(j); }

		if( B_SE(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_SE at right boundary in FG"));  }

		if( B_SW(flags.at(imax+1).at(j))){ F.at(imax).at(j) = U.at(imax).at(j); G.at(imax+1).at(j-1) = V.at(imax+1).at(j-1); }

		// Inflow - only x direction
		if((flags.at(imax+1).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at right boundary in FG")); }

		// Outflow - only x direction
		if((flags.at(imax+1).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(j) = U.at(imax).at(j);}

	}

	// ------------------------------------------------------- Top boundary ------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		if( B_N(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top boundary in FG")); }

		if( B_S(flags.at(i).at(jmax+1)) ){ G.at(i).at(jmax) = V.at(i).at(jmax); }

		if( B_W(flags.at(i).at(jmax+1)) ){ F.at(i-1).at(jmax+1) = U.at(i-1).at(jmax+1); }

		if( B_E(flags.at(i).at(jmax+1)) ){ F.at(i).at(jmax+1) = U.at(i).at(jmax+1); }

		if( B_NE(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NE at top boundary in FG")); }

		if( B_NW(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top boundary in FG")); }

		if( B_SE(flags.at(i).at(jmax+1)) ){ F.at(i).at(jmax+1) = U.at(i).at(jmax+1); G.at(i).at(jmax) = V.at(i).at(jmax); }

		if( B_SW(flags.at(i).at(jmax+1))){ F.at(i-1).at(jmax+1) = U.at(i-1).at(jmax+1); G.at(i).at(jmax) = V.at(i).at(jmax); }

		// Inflow - only x direction
		if((flags.at(i).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(jmax+1) = U.at(i).at(jmax+1);}

		// Outflow - only x direction
		if((flags.at(i).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(jmax+1) = U.at(i-1).at(jmax+1);}

	}

	// ---------------------------------------------------- Bottom boundary ------------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		if( B_N(flags.at(i).at(0)) ){ G.at(i).at(0) = V.at(i).at(0); }

		if( B_S(flags.at(i).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom boundary in FG")); }

		if( B_W(flags.at(i).at(0)) ){ F.at(i-1).at(0) = U.at(i-1).at(0); }

		if( B_E(flags.at(i).at(0)) ){ F.at(i).at(0) = U.at(i).at(0); }

		if( B_NE(flags.at(i).at(0)) ){ F.at(i).at(0) = U.at(i).at(0); G.at(i).at(0) = V.at(i).at(0); }

		if( B_NW(flags.at(i).at(0)) ){ F.at(i-1).at(0) = U.at(i-1).at(0); G.at(i).at(0) = V.at(i).at(0); }

		if( B_SE(flags.at(i).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary in FG")); }

		if( B_SW(flags.at(i).at(0))){ throw std::runtime_error(std::string("impossible to have B_WW at bottom boundary in FG")); }

		// Inflow - only x direction
		if((flags.at(i).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(0) = U.at(i).at(0);}

		// Outflow - only x direction
		if((flags.at(i).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(0) = U.at(i-1).at(0);}

	}

	// --------------------------------------------------- Top-left corner --------------------------------------------------------
	if( B_N(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-left boundary in FG")); }

	if( B_S(flags.at(0).at(jmax+1)) ){ G.at(0).at(jmax) = V.at(0).at(jmax); }

	if( B_W(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_W at top-left boundary in FG")); }

	if( B_E(flags.at(0).at(jmax+1)) ){ F.at(0).at(jmax+1) = U.at(0).at(jmax+1); }

	if( B_NE(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary in FG")); }

	if( B_NW(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary in FG")); }

	if( B_SE(flags.at(0).at(jmax+1)) ){ F.at(0).at(jmax+1) = U.at(0).at(jmax+1); G.at(0).at(jmax) = V.at(0).at(jmax); }

	if( B_SW(flags.at(0).at(jmax+1))){ throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(0).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(jmax+1) = U.at(0).at(jmax+1);}

	// Outflow - only x direction
	if((flags.at(0).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at top-left boundary in FG"));}


	// ---------------------------------------------------- Top-right corner -----------------------------------------------------------
	if( B_N(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-right boundary in FG")); }

	if( B_S(flags.at(imax+1).at(jmax+1)) ){ G.at(imax+1).at(jmax) = V.at(imax+1).at(jmax); }

	if( B_W(flags.at(imax+1).at(jmax+1)) ){ F.at(imax).at(jmax+1) = U.at(imax).at(jmax+1); }

	if( B_E(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_W at top-right boundary in FG")); }

	if( B_NE(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-right boundary in FG")); }

	if( B_NW(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary in FG")); }

	if( B_SE(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary in FG")); }

	if( B_SW(flags.at(imax+1).at(jmax+1))){ F.at(imax).at(jmax+1) = U.at(imax).at(jmax+1); G.at(imax+1).at(jmax) = V.at(imax+1).at(jmax); }

	// Inflow - only x direction
	if((flags.at(imax+1).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at top-right boundary in FG"));}

	// Outflow - only x direction
	if((flags.at(imax+1).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(jmax+1) = U.at(imax).at(jmax+1);}


	// ----------------------------------------------------- Bottom-left corner -----------------------------------------------------------

	if( B_N(flags.at(0).at(0)) ){ G.at(0).at(0) = V.at(0).at(0); }

	if( B_S(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary in FG")); }

	if( B_W(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary in FG")); }

	if( B_E(flags.at(0).at(0)) ){ F.at(0).at(0) = U.at(0).at(0); }

	if( B_NE(flags.at(0).at(0)) ){ F.at(0).at(0) = U.at(0).at(0); G.at(0).at(0) = V.at(0).at(0); }

	if( B_NW(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary in FG")); }

	if( B_SE(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary in FG")); }

	if( B_SW(flags.at(0).at(0))){ throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(0).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(0) = U.at(0).at(0);}

	// Outflow - only x direction
	if((flags.at(0).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at bottom-left boundary in FG"));}

	// ------------------------------------------------------ Bottom-right corner ------------------------------------------------------------

	if( B_N(flags.at(imax+1).at(0)) ){ G.at(imax+1).at(0) = V.at(imax+1).at(0); }

	if( B_S(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary in FG")); }

	if( B_W(flags.at(imax+1).at(0)) ){ F.at(imax).at(0) = U.at(imax).at(0); }

	if( B_E(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary in FG")); }

	if( B_NE(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary in FG")); }

	if( B_NW(flags.at(imax+1).at(0)) ){ F.at(imax).at(0) = U.at(imax).at(0); G.at(imax+1).at(0) = V.at(imax+1).at(0); }

	if( B_SE(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary in FG")); }

	if( B_SW(flags.at(imax+1).at(0))){ throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(imax+1).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at bottom-right boundary in FG"));}

	// Outflow - only x direction
	if((flags.at(imax+1).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(0) = U.at(imax).at(0);}


	// ----------------------------------------- Inner cells -------------------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		for(size_t j = 1; j <= jmax; j++){


			if( B_N(flags.at(i).at(j)) ){ G.at(i).at(j) = V.at(i).at(j); }

			if( B_S(flags.at(i).at(j)) ){ G.at(i).at(j-1) = V.at(i).at(j-1); }

			if( B_W(flags.at(i).at(j)) ){ F.at(i-1).at(j) = U.at(i-1).at(j); }

			if( B_E(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); }

			if( B_NE(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); G.at(i).at(j) = V.at(i).at(j); }

			if( B_NW(flags.at(i).at(j)) ){ F.at(i-1).at(j) = U.at(i-1).at(j); G.at(i).at(j) = V.at(i).at(j); }

			if( B_SE(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); G.at(i).at(j-1) = V.at(i).at(j-1); }

			if( B_SW(flags.at(i).at(j))){ F.at(i-1).at(j) = U.at(i-1).at(j); G.at(i).at(j-1) = V.at(i).at(j-1); }

			// Inflow - only x direction
			if((flags.at(i).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(j) = U.at(i).at(j);}

			// Outflow - only x direction
			if((flags.at(i).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(j) = U.at(i-1).at(j);}

		}
	}





	// ------------------------------------------------------ F and G calculations -------------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		for( size_t j = 1; j <= jmax; j++){

			if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE)){

				F.at(i).at(j) =  (

						(1.0 / Re) *( ( ( U.at(i-1).at(j) - 2.0 * U.at(i).at(j) + U.at(i+1).at(j) ) / ( dx * dx ) ) +  ( ( U.at(i).at(j-1) - 2.0 * U.at(i).at(j) + U.at(i).at(j+1) ) / ( dy * dy ) ) )

						// d(u^2)/d(x)
						- (0.25) * ( 1.0 / dx ) * (  pow( ( U.at(i).at(j) + U.at(i+1).at(j) ), 2.0 )  -   pow( ( U.at(i-1).at(j) + U.at(i).at(j) ), 2.0 )
													 + alpha * ( ( fabs( U.at(i).at(j) + U.at(i+1).at(j) )  * ( U.at(i).at(j) - U.at(i+1).at(j)) ) -   (  fabs( U.at(i-1).at(j) + U.at(i).at(j)) * ( U.at(i-1).at(j) - U.at(i).at(j) ) ) )
						)

						// d(uv)/d(y)
						- (0.25) * ( 1.0 / dy ) * ( ( V.at(i).at(j) + V.at(i+1).at(j) ) * ( U.at(i).at(j) + U.at(i).at(j+1) ) - ( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) + U.at(i).at(j) )
													+ alpha * ( fabs(V.at(i).at(j) + V.at(i+1).at(j)) * ( U.at(i).at(j) - U.at(i).at(j+1) ) - fabs( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) - U.at(i).at(j) )   )
						)

						+ GX);

			}
		}
	}

	for(size_t i = 1; i <= imax; i++){

		for( size_t j = 1; j <= jmax; j++){

			if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE))
			{
				G.at(i).at(j) =  (
						(1.0 / Re) * (((V.at(i - 1).at(j) - 2.0 * V.at(i).at(j) + V.at(i + 1).at(j)) / (dx * dx)) +
									  ((V.at(i).at(j - 1) - 2.0 * V.at(i).at(j) + V.at(i).at(j + 1)) / (dy * dy)))

						// d(v^2) / dy
						- (0.25) * (1.0 / dy) *
						  (pow((V.at(i).at(j) + V.at(i).at(j + 1)), 2.0) - pow((V.at(i).at(j - 1) + V.at(i).at(j)), 2.0)
						   + alpha * ((fabs(V.at(i).at(j) + V.at(i).at(j + 1)) * (V.at(i).at(j) - V.at(i).at(j + 1))) -
									  (fabs(V.at(i).at(j - 1) + V.at(i).at(j)) * (V.at(i).at(j - 1) - V.at(i).at(j))))
						  )

						// d(uv)/dx
						- (0.25) * (1.0 / dx) *
						  ((U.at(i).at(j) + U.at(i).at(j + 1)) * (V.at(i).at(j) + V.at(i + 1).at(j)) -
						   (U.at(i - 1).at(j) + U.at(i - 1).at(j + 1)) * (V.at(i - 1).at(j) + V.at(i).at(j))
						   + alpha * (fabs(U.at(i).at(j) + U.at(i).at(j + 1)) * (V.at(i).at(j) - V.at(i + 1).at(j)) -
									  fabs(U.at(i - 1).at(j) + U.at(i - 1).at(j + 1)) *
									  (V.at(i - 1).at(j) - V.at(i).at(j)))
						  )

						+ GY);

			}
		}
	}





}

// This operation computes the right hand side of the pressure poisson equation.
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
        matrix<double> &RS)
{

    for(size_t i = 1; i <= imax ; i++){

        for( size_t j = 1; j <= jmax ; j++)

        {

            RS.at(i).at(j) = (1.0 / dt) *
                             ((U_star.at(i).at(j) - U_star.at(i-1).at(j)) / dx + (V_star.at(i).at(j) - V_star.at(i).at(j-1)) / dy) +
                             (F.at(i).at(j) - F.at(i-1).at(j)) / dx + (G.at(i).at(j) - G.at(i).at(j-1)) / dy;

        }

    }

}

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
){

	for(size_t i = 1; i <= imax ; i++){

		for( size_t j = 1; j <= jmax ; j++)
		{
			if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE))
			{
				RS.at(i).at(j) = (1.0 / dt) *
				                             ((U_star.at(i).at(j) - U_star.at(i-1).at(j)) / dx + (V_star.at(i).at(j) - V_star.at(i).at(j-1)) / dy) +
				                             (F.at(i).at(j) - F.at(i-1).at(j)) / dx + (G.at(i).at(j) - G.at(i).at(j-1)) / dy;
			}
		}
	}


}


// Determines the maximal time step size
void calculate_dt(double Re, double tau, double *dt, double dx, double dy, int imax, int jmax, const matrix<double> &U, const matrix<double> &V)
{

        double u_max = maxElementOfMatrix(U);
        double v_max = maxElementOfMatrix(V);

        double term1 = (Re / 2.0) * ( 1.0 / ((1.0 / (dx * dx) ) + ( 1.0 / (dy * dy) ) ) );
        double term2 = dx / fabs(u_max);
        double term3 = dy / fabs(v_max);

        double min1 = std::min(term1,term2);
        double min2 = std::min(min1,term3);


        //*dt = tau * min2;
        if( (tau > 0) && (tau < 1))
         {
           *dt = tau * min2;
         }

}

void calculate_dt_arbitrary(double Re,
        double tau,
        double *dt,
        double dx,
        double dy,
        int imax,
        int jmax,
		const matrix<double> &U,
		const matrix<double> &V){


	double u_max = maxElementOfMatrix(U);
	double v_max = maxElementOfMatrix(V);

	double term = ( 1.0 / ((1.0 / (dx * dx) ) + ( 1.0 / (dy * dy) ) ) ) ;
	double term1 = (Re / 2.0) * term;
	double term2 = dx / fabs(u_max);
	double term3 = dy / fabs(v_max);


	double min1 = std::min(term1,term2);
	double min2 = std::min(min1,term3);



	//*dt = tau * min2;
	if( (tau > 0) && (tau < 1))
	{
		*dt = tau * min2;
	}

}

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
        )
{


    if(*method == "ExplicitEuler"){

        for( size_t i = 1; i <= imax-1; i++){

            for( size_t j = 1; j <= jmax ; j++){

                U.at(i).at(j) = U.at(i).at(j) + dt * UU.at(0).at(i).at(j);

            }

        }

        for( size_t i = 1; i <= imax; i++){

            for( size_t j = 1; j <= jmax-1 ; j++) {

                V.at(i).at(j) =  V.at(i).at(j) + dt * VV.at(0).at(i).at(j);

            }

        }

    }

    else if(*method == "Heun"){

        for( size_t i = 1; i <= imax-1; i++){

            for( size_t j = 1; j <= jmax ; j++){

                U.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * (UU.at(0).at(i).at(j) + UU.at(1).at(i).at(j) );

            }

        }

        for( size_t i = 1; i <= imax; i++){

            for( size_t j = 1; j <= jmax-1 ; j++) {

                V.at(i).at(j) =  V.at(i).at(j) + dt * 0.5 * (VV.at(0).at(i).at(j) + VV.at(1).at(i).at(j) );

            }

        }

    }
    else if(*method == "RungeKutta")
	{
		for( size_t i = 1; i <= imax-1; i++){

			for( size_t j = 1; j <= jmax ; j++){

				U.at(i).at(j) = U.at(i).at(j) + dt * (  (UU.at(0).at(i).at(j) + 2.0 * UU.at(1).at(i).at(j) + 2.0 * UU.at(2).at(i).at(j) + UU.at(3).at(i).at(j) ) / 6.0);

			}

		}


		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax-1 ; j++){

				V.at(i).at(j) = V.at(i).at(j) + dt * (  (VV.at(0).at(i).at(j) + 2.0 * VV.at(1).at(i).at(j) + 2.0 * VV.at(2).at(i).at(j) + VV.at(3).at(i).at(j) ) / 6.0);

			}

		}
	}

}

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

){

	if(*method == "ExplicitEuler"){

		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++){

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				U.at(i).at(j) = U.at(i).at(j) + dt * UU.at(0).at(i).at(j);
				}
			}

		}

		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++) {

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				V.at(i).at(j) =  V.at(i).at(j) + dt * VV.at(0).at(i).at(j);
				}

			}

		}

	}

	else if(*method == "Heun"){

		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++){

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				U.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * (UU.at(0).at(i).at(j) + UU.at(1).at(i).at(j) );
				}

			}

		}

		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++) {

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				V.at(i).at(j) =  V.at(i).at(j) + dt * 0.5 * (VV.at(0).at(i).at(j) + VV.at(1).at(i).at(j) );
				}

			}

		}

	}
	else if(*method == "RungeKutta")
	{
		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++){

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				U.at(i).at(j) = U.at(i).at(j) + dt * (  (UU.at(0).at(i).at(j) + 2.0 * UU.at(1).at(i).at(j) + 2.0 * UU.at(2).at(i).at(j) + UU.at(3).at(i).at(j) ) / 6.0);
				}

			}

		}


		for( size_t i = 1; i <= imax; i++){

			for( size_t j = 1; j <= jmax ; j++){

				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				V.at(i).at(j) = V.at(i).at(j) + dt * (  (VV.at(0).at(i).at(j) + 2.0 * VV.at(1).at(i).at(j) + 2.0 * VV.at(2).at(i).at(j) + VV.at(3).at(i).at(j) ) / 6.0);
				}

			}

		}
	}






}



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
){


    if (method == "ExplicitEuler")
    {
        for(size_t i = 0; i <= imax+1; i++)
        {

            for(size_t j = 0; j <= jmax+1; j++)
            {

                U_star.at(i).at(j) = U.at(i).at(j);

            }

        }


        for(size_t i = 0; i <= imax+1; i++)
        {

            for(size_t j = 0; j <= jmax+1; j++)
            {

                V_star.at(i).at(j) = V.at(i).at(j);


            }

        }




    }
    else if(method == "Heun")
    {


        if(k == 0)
        {


            for(size_t i = 0; i <= imax+1; i++)
            {

                for(size_t j = 0; j <= jmax+1; j++)
                {

                    U_star.at(i).at(j) = U.at(i).at(j);

                }

            }


            for(size_t i = 0; i <= imax+1; i++)
            {

                for(size_t j = 0; j <= jmax+1; j++)
                {

                    V_star.at(i).at(j) = V.at(i).at(j);


                }

            }

        }

        else if( k == 1){

            for(size_t i = 0; i <= imax+1; i++)

            {

                for(size_t j = 0; j <= jmax+1; j++)

                {

                    U_star.at(i).at(j) = U.at(i).at(j) + dt * UU.at(0).at(i).at(j);

                }

            }

            for(size_t i = 0; i <= imax+1; i++)

            {

                for(size_t j = 0; j <= jmax+1; j++)

                {

                    V_star.at(i).at(j) = V.at(i).at(j) + dt * VV.at(0).at(i).at(j);

                }

            }


        }


    }
    else if( method == "RungeKutta")
       {


           if(k == 0)
           {


               for(size_t i = 0; i <= imax+1; i++)
               {

                   for(size_t j = 0; j <= jmax+1; j++)
                   {

                       U_star.at(i).at(j) = U.at(i).at(j);

                   }

               }


               for(size_t i = 0; i <= imax+1; i++)
               {

                   for(size_t j = 0; j <= jmax+1; j++)
                   {

                       V_star.at(i).at(j) = V.at(i).at(j);


                   }

               }

           }
           else if( k == 1){

               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       U_star.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * UU.at(0).at(i).at(j);

                   }

               }

               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       V_star.at(i).at(j) = V.at(i).at(j) + dt * 0.5 * VV.at(0).at(i).at(j);

                   }

               }


           }
           else if( k == 2)
           {

               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       U_star.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * UU.at(1).at(i).at(j);

                   }

               }

               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       V_star.at(i).at(j) = V.at(i).at(j) + dt * 0.5 * VV.at(1).at(i).at(j);

                   }

               }

           }
           else if(k == 3)
           {
               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       U_star.at(i).at(j) = U.at(i).at(j) + dt * UU.at(2).at(i).at(j);

                   }

               }

               for(size_t i = 0; i <= imax+1; i++)

               {

                   for(size_t j = 0; j <= jmax+1; j++)

                   {

                       V_star.at(i).at(j) = V.at(i).at(j) + dt * VV.at(2).at(i).at(j);

                   }

               }


           }






       }


}

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
){

	if (method == "ExplicitEuler")
	{
		for(size_t i = 0; i <= imax+1; i++)
		{

			for(size_t j = 0; j <= jmax+1; j++)
			{
				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				U_star.at(i).at(j) = U.at(i).at(j);
				}
			}

		}


		for(size_t i = 0; i <= imax+1; i++)
		{

			for(size_t j = 0; j <= jmax+1; j++)
			{
				if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
				V_star.at(i).at(j) = V.at(i).at(j);
				}

			}

		}




	}
	else if(method == "Heun")
	{


		if(k == 0)
		{


			for(size_t i = 0; i <= imax+1; i++)
			{

				for(size_t j = 0; j <= jmax+1; j++)
				{

					if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					U_star.at(i).at(j) = U.at(i).at(j);
					}
				}

			}


			for(size_t i = 0; i <= imax+1; i++)
			{

				for(size_t j = 0; j <= jmax+1; j++)
				{
					if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					V_star.at(i).at(j) = V.at(i).at(j);
					}

				}

			}

		}

		else if( k == 1){

			for(size_t i = 0; i <= imax+1; i++)

			{

				for(size_t j = 0; j <= jmax+1; j++)

				{
					if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					U_star.at(i).at(j) = U.at(i).at(j) + dt * UU.at(0).at(i).at(j);
					}
				}

			}

			for(size_t i = 0; i <= imax+1; i++)

			{

				for(size_t j = 0; j <= jmax+1; j++)

				{
					if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					V_star.at(i).at(j) = V.at(i).at(j) + dt * VV.at(0).at(i).at(j);
					}
				}

			}


		}


	}
	else if( method == "RungeKutta")
	   {


		   if(k == 0)
		   {


			   for(size_t i = 0; i <= imax+1; i++)
			   {

				   for(size_t j = 0; j <= jmax+1; j++)
				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   U_star.at(i).at(j) = U.at(i).at(j);
					   }
				   }

			   }


			   for(size_t i = 0; i <= imax+1; i++)
			   {

				   for(size_t j = 0; j <= jmax+1; j++)
				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   V_star.at(i).at(j) = V.at(i).at(j);
					   }

				   }

			   }

		   }
		   else if( k == 1){

			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   U_star.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * UU.at(0).at(i).at(j);
					   }
				   }

			   }

			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   V_star.at(i).at(j) = V.at(i).at(j) + dt * 0.5 * VV.at(0).at(i).at(j);
					   }
				   }

			   }


		   }
		   else if( k == 2)
		   {

			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   U_star.at(i).at(j) = U.at(i).at(j) + dt * 0.5 * UU.at(1).at(i).at(j);
					   }
				   }

			   }

			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   V_star.at(i).at(j) = V.at(i).at(j) + dt * 0.5 * VV.at(1).at(i).at(j);
					   }
				   }

			   }

		   }
		   else if(k == 3)
		   {
			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   U_star.at(i).at(j) = U.at(i).at(j) + dt * UU.at(2).at(i).at(j);
					   }
				   }

			   }

			   for(size_t i = 0; i <= imax+1; i++)

			   {

				   for(size_t j = 0; j <= jmax+1; j++)

				   {
					   if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
					   V_star.at(i).at(j) = V.at(i).at(j) + dt * VV.at(2).at(i).at(j);
					   }
				   }

			   }


		   }


	   }



}



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
){


    for (size_t i = 1; i <= imax-1; ++i)

    {

        for (size_t j = 1; j <= jmax; ++j)

        {

            UU.at(k).at(i).at(j) = F.at(i).at(j) - 1.0 / dx * (P.at(i+1).at(j) - P.at(i).at(j));


        }

    }


    for (size_t i = 1; i <= imax; ++i)

    {
        for (size_t j = 1; j <= jmax-1; ++j)

        {

            VV.at(k).at(i).at(j) = G.at(i).at(j) - 1.0 / dy * (P.at(i).at(j + 1) - P.at(i).at(j));

        }
    }

}


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
){

	for (size_t i = 1; i <= imax-1; ++i)

	{

		for (size_t j = 1; j <= jmax; ++j)

		{
			if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
			UU.at(k).at(i).at(j) = F.at(i).at(j) - 1.0 / dx * (P.at(i+1).at(j) - P.at(i).at(j));
			}

		}

	}


	for (size_t i = 1; i <= imax; ++i)

	{
		for (size_t j = 1; j <= jmax-1; ++j)

		{
			if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
			VV.at(k).at(i).at(j) = G.at(i).at(j) - 1.0 / dy * (P.at(i).at(j + 1) - P.at(i).at(j));
			}
		}
	}

}

