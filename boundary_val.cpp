#include "boundary_val.hpp"
#include "datastructures.hpp"

void boundaryvalues(int imax, int jmax, matrix<double> &U, matrix<double> &V ) {

    for( size_t j = 1; j  <= jmax; j++){
    	// U component at Vertical walls = 0
        U.at(0).at(j) = 0; // Set cell
        U.at(imax).at(j) = 0;

        // V component at Vertical walls
        V.at(0).at(j) = -V.at(1).at(j);
        V.at(imax+1).at(j) = -V.at(imax).at(j);

    }

    for( size_t i = 1; i <= imax; i++){

    	// V component at Horizontal walls = 0
    	V.at(i).at(0) = 0;
    	V.at(i).at(jmax) = 0;

    	// U component at Bottom wall
    	U.at(i).at(0) = - U.at(i).at(1);
    	// U component at Top wall
    	U.at(i).at(jmax+1) = 2-U.at(i).at(jmax);
    }



}

void boundaryvalues_arbitrary(int imax, int jmax, matrix<double> &U, matrix<double> &V, const matrix<unsigned int> &flag){

	// ---------------------------------------------------- left boundary --------------------------------------------------------------
	for(int j = 1; j<=jmax; ++j){

		switch(flag.at(0).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW )){
			case flag_10bit::NO_SLIP: //No slip conditions
				if (B_N(flag.at(0).at(j)) == 1){
					V.at(0).at(j) = 0;
					U.at(0).at(j) = -U.at(0).at(j+1);
				}

				if (B_W(flag.at(0).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_W at left boundary "));
				}

				if (B_S(flag.at(0).at(j)) == 1){
					V.at(0).at(j-1) = 0;
					U.at(0).at(j) = -U.at(0).at(j-1);
				}

				if (B_E(flag.at(0).at(j)) == 1){
					U.at(0).at(j) = 0;
					V.at(0).at(j-1) = -V.at(1).at(j-1);
					V.at(0).at(j) = -V.at(1).at(j);
				}

				if (B_NE(flag.at(0).at(j)) == 1){
					U.at(0).at(j) = 0;
					V.at(0).at(j) = 0;
					V.at(0).at(j-1) = -V.at(1).at(j-1);
				}

				if (B_NW(flag.at(0).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
				}

				if (B_SE(flag.at(0).at(j)) == 1){
					U.at(0).at(j)=0;
					V.at(0).at(j-1)=0;
					V.at(0).at(j) = -V.at(1).at(j);
				}

				if (B_SW(flag.at(0).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
				}
				break;
			case flag_10bit::FREE_SLIP://Free Slip conditions
				if (B_E(flag.at(0).at(j)) == 1){
					U.at(0).at(j) = 0;
					V.at(0).at(j) = V.at(1).at(j);
					V.at(0).at(j-1) = V.at(1).at(j-1);
				}

				if (B_W(flag.at(0).at(j)) == 1) {
					throw std::runtime_error(std::string("impossible to have B_W at left boundary "));
				}

				if (B_N(flag.at(0).at(j)) == 1) {
					V.at(0).at(j) = 0;
					U.at(0).at(j) = U.at(0).at(j+1);
				}

				if (B_S(flag.at(0).at(j)) == 1){
					V.at(0).at(j-1) = 0;
					U.at(0).at(j) = U.at(0).at(j-1);
				}

				if (B_NE(flag.at(0).at(j)) == 1){
					U.at(0).at(j) = 0;
					V.at(0).at(j) = 0;
					V.at(0).at(j-1) = V.at(1).at(j-1);
				}

				if (B_NW(flag.at(0).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
				}

				if (B_SE(flag.at(0).at(j)) == 1){
					U.at(0).at(j)=0;
					V.at(0).at(j-1)=0;
					V.at(0).at(j) = V.at(1).at(j);
				}

				if (B_SW(flag.at(0).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SW at left boundary "));
				}
				break;

			case flag_10bit::OUTFLOW: // assuming only x-direction
				throw std::runtime_error(std::string("impossible to have outflow at left boundary (only x direction assumed"));
				break;

			case flag_10bit::INFLOW: // assuming only x-direction

				U.at(0).at(j)=1;
				V.at(0).at(j) = 0;
				V.at(0).at(j-1) = 0;
				break;

			default:
				break;

		}

	}

	// ------------------------------------------------------- Right Boundary ----------------------------------------------------------------
	for(int j = 1; j<=jmax; ++j){

		switch(flag.at(imax+1).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
			case flag_10bit::NO_SLIP: //No slip conditions
				if (B_N(flag.at(imax+1).at(j)) == 1){
					V.at(imax+1).at(j) = 0;
					U.at(imax).at(j) = -U.at(imax).at(j+1);
					U.at(imax+1).at(j) = -U.at(imax+1).at(j+1);
				}

				if (B_W(flag.at(imax+1).at(j)) == 1){
					U.at(imax).at(j) = 0;
					V.at(imax+1).at(j-1) = -V.at(imax).at(j-1);
					V.at(imax+1).at(j) = -V.at(imax).at(j);
				}

				if (B_S(flag.at(imax+1).at(j)) == 1){
					V.at(imax+1).at(j-1) = 0;
					U.at(imax).at(j) = -U.at(imax).at(j-1);
					U.at(imax+1).at(j) = -U.at(imax+1).at(j-1);
				}

				if (B_E(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_E at right boundary "));
				}

				if (B_NE(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NE at right boundary "));
				}

				if (B_NW(flag.at(imax+1).at(j)) == 1){
					U.at(imax).at(j) = 0;
					U.at(imax+1).at(j) = - U.at(imax+1).at(j+1);
					V.at(imax+1).at(j) = 0;
					V.at(imax+1).at(j-1) = -V.at(imax).at(j-1);
				}

				if (B_SE(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SE at right boundary "));
				}

				if (B_SW(flag.at(imax+1).at(j)) == 1){
					U.at(imax).at(j) = 0;
					U.at(imax+1).at(j) = -U.at(imax+1).at(j-1);
					V.at(imax+1).at(j-1) = 0;
					V.at(imax+1).at(j) = -V.at(imax).at(j);
				}
				break;

			case flag_10bit::FREE_SLIP://Free Slip conditions
				if (B_E(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_E at right boundary "));
				}

				if (B_W(flag.at(imax+1).at(j)) == 1) {
					U.at(imax).at(j) = 0;
					V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
					V.at(imax+1).at(j) = V.at(imax).at(j);
				}

				if (B_N(flag.at(imax+1).at(j)) == 1) {
					V.at(imax+1).at(j) = 0;
					U.at(imax).at(j) = U.at(imax).at(j+1);
					U.at(imax+1).at(j) = U.at(imax+1).at(j+1);
				}

				if (B_S(flag.at(imax+1).at(j)) == 1){
					V.at(imax+1).at(j-1) = 0;
					U.at(imax).at(j) = U.at(imax).at(j-1);
					U.at(imax+1).at(j) = U.at(imax+1).at(j-1);
				}

				if (B_NE(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NE at right boundary "));
				}

				if (B_NW(flag.at(imax+1).at(j)) == 1){
					U.at(imax).at(j) = 0;
					U.at(imax+1).at(j) = U.at(imax+1).at(j+1);
					V.at(imax+1).at(j) = 0;
					V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
				}

				if (B_SE(flag.at(imax+1).at(j)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SE at right boundary "));
				}

				if (B_SW(flag.at(imax+1).at(j)) == 1){
					U.at(imax).at(j) = 0;
					U.at(imax+1).at(j) = U.at(imax+1).at(j-1);
					V.at(imax+1).at(j-1) = 0;
					V.at(imax+1).at(j) = V.at(imax).at(j);
				}
				break;

			case flag_10bit::OUTFLOW: // assuming only x-direction
				U.at(imax+1).at(j) = U.at(imax).at(j);
				V.at(imax+1).at(j) = V.at(imax).at(j);
				V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
				break;

			case flag_10bit::INFLOW: // assuming only x-direction
				U.at(imax+1).at(j)=1;
				V.at(imax+1).at(j) = 0;
				V.at(imax+1).at(j-1) = 0;
				break;

			default:
				break;
		}
	}

	//  ---------------------------------------------------- Top Boundary ----------------------------------------------------------------
	for(int i = 1; i<=imax; ++i){
		switch(flag.at(i).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
			case flag_10bit::NO_SLIP: //No slip conditions
				if (B_N(flag.at(i).at(jmax+1)) == 1){
					throw std::runtime_error(std::string("impossible to have B_N at top boundary "));
				}

				if (B_W(flag.at(i).at(jmax+1)) == 1){
					U.at(i-1).at(jmax+1) = 0;
					V.at(i).at(jmax) = -V.at(i-1).at(jmax);
					V.at(i).at(jmax+1) = -V.at(i-1).at(jmax+1);
				}

				if (B_S(flag.at(i).at(jmax+1)) == 1){
					V.at(i).at(jmax) = 0;
					U.at(i-1).at(jmax+1) = -U.at(i-1).at(jmax);
					U.at(i).at(jmax+1) = -U.at(i).at(jmax);
				}

				if (B_E(flag.at(i).at(jmax+1)) == 1){
					U.at(i).at(jmax+1) = 0;
					V.at(i).at(jmax) = -V.at(i+1).at(jmax);
					V.at(i).at(jmax+1) = -V.at(i+1).at(jmax+1);
				}

				if (B_NE(flag.at(i).at(jmax+1)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NE at top boundary "));
				}

				if (B_NW(flag.at(i).at(jmax+1)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NW at top boundary "));
				}

				if (B_SE(flag.at(i).at(jmax+1)) == 1){
					U.at(i).at(jmax+1)=0;
					U.at(i-1).at(jmax+1) = -U.at(i-1).at(jmax);
					V.at(i).at(jmax)=0;
					V.at(i).at(jmax+1) = -V.at(i+1).at(jmax+1);
				}

				if (B_SW(flag.at(i).at(jmax+1)) == 1){
					U.at(i-1).at(jmax+1) = 0;
					U.at(i).at(jmax+1) = -U.at(i).at(jmax);
					V.at(i).at(jmax) = 0;
					V.at(i).at(jmax+1) = -V.at(i-1).at(jmax+1);
				}
				break;

			case flag_10bit::FREE_SLIP://Free Slip conditions
				if (B_E(flag.at(i).at(jmax+1)) == 1){
					U.at(i).at(jmax+1) = 0;
					V.at(i).at(jmax+1) = V.at(i+1).at(jmax+1);
					V.at(i).at(jmax) = V.at(i+1).at(jmax);
				}

				if (B_W(flag.at(i).at(jmax+1)) == 1) {
					U.at(i-1).at(jmax+1) = 0;
					V.at(i).at(jmax) = V.at(i-1).at(jmax);
					V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
				}

				if (B_N(flag.at(i).at(jmax+1)) == 1) {
					throw std::runtime_error(std::string("impossible to have B_N at top boundary "));
				}

				if (B_S(flag.at(i).at(jmax+1)) == 1){
					V.at(i).at(jmax) = 0;
					U.at(i-1).at(jmax+1) = U.at(i-1).at(jmax);
					U.at(i).at(jmax+1) = U.at(i).at(jmax);
				}

				if (B_NE(flag.at(i).at(jmax+1)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NE at top boundary "));
				}

				if (B_NW(flag.at(i).at(jmax+1)) == 1){
					throw std::runtime_error(std::string("impossible to have B_NW at top boundary "));
				}

				if (B_SE(flag.at(i).at(jmax+1)) == 1){
					U.at(i).at(jmax+1)=0;
					U.at(i-1).at(jmax+1) = U.at(i-1).at(jmax);
					V.at(i).at(jmax)=0;
					V.at(i).at(jmax+1) = V.at(i+1).at(jmax+1);
				}

				if (B_SW(flag.at(i).at(jmax+1)) == 1){
					U.at(i-1).at(jmax+1) = 0;
					U.at(i).at(jmax+1) = U.at(i).at(jmax);
					V.at(i).at(jmax) = 0;
					V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
				}
				break;

			case flag_10bit::OUTFLOW: // assuming only x-direction
				U.at(i).at(jmax+1) = U.at(i-1).at(jmax+1);
				V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
				V.at(i).at(jmax) = V.at(i-1).at(jmax);
				break;

			case flag_10bit::INFLOW: // assuming only x-direction
				U.at(i).at(jmax+1)=1;
				V.at(i).at(jmax+1) = 0;
				V.at(i).at(jmax) = 0;
				break;

			default:
				break;
		}

	}

	// ------------------------------------------------- Bottom Boundary -----------------------------------------------------------------
	for(int i = 1; i<=imax; ++i){
		switch(flag.at(i).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
			case flag_10bit::NO_SLIP: //No slip conditions
				if (B_N(flag.at(i).at(0)) == 1){
					V.at(i).at(0) = 0;
					U.at(i-1).at(0) = -U.at(i-1).at(1);
					U.at(i).at(0) = -U.at(i).at(1);
				}

				if (B_W(flag.at(i).at(0)) == 1){
					U.at(i-1).at(0) = 0;
					V.at(i).at(0) = -V.at(i-1).at(0);
				}

				if (B_S(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_S at bottom boundary "));
				}

				if (B_E(flag.at(i).at(0)) == 1){
					U.at(i).at(0) = 0;
					V.at(i).at(0) = -V.at(i+1).at(0);
				}

				if (B_NE(flag.at(i).at(0)) == 1){
					U.at(i).at(0) = 0;
					U.at(i-1).at(0) = -U.at(i-1).at(1);
					V.at(i).at(0) = 0;
				}

				if (B_NW(flag.at(i).at(0)) == 1){
					U.at(i-1).at(0) = 0;
					U.at(i).at(0) = - U.at(i).at(1);
					V.at(i).at(0) = 0;
				}

				if (B_SE(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary "));
				}

				if (B_SW(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SW at bottom boundary "));
				}
				break;

			case flag_10bit::FREE_SLIP://Free Slip conditions
				if (B_E(flag.at(i).at(0)) == 1){
					U.at(i).at(0) = 0;
					V.at(i).at(0) = V.at(i+1).at(0);
				}

				if (B_W(flag.at(i).at(0)) == 1) {
					U.at(i-1).at(0) = 0;
					V.at(i).at(0) = V.at(i-1).at(0);
				}

				if (B_N(flag.at(i).at(0)) == 1) {
					V.at(i).at(0) = 0;
					U.at(i-1).at(0) = U.at(i-1).at(1);
					U.at(i).at(0) = U.at(i).at(1);
				}

				if (B_S(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_S at top boundary "));
				}

				if (B_NE(flag.at(i).at(0)) == 1){
					U.at(i).at(0) = 0;
					U.at(i-1).at(0) = U.at(i-1).at(1);
					V.at(i).at(0) = 0;
				}

				if (B_NW(flag.at(i).at(0)) == 1){
					U.at(i-1).at(0) = 0;
					U.at(i).at(0) = U.at(i).at(1);
					V.at(i).at(0) = 0;
				}

				if (B_SE(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SE at top boundary "));
				}

				if (B_SW(flag.at(i).at(0)) == 1){
					throw std::runtime_error(std::string("impossible to have B_SW at top boundary "));
				}
				break;

			case flag_10bit::OUTFLOW: // assuming only x-direction
				U.at(i).at(0) = U.at(i-1).at(0);
				V.at(i).at(0) = V.at(i-1).at(0);
				break;

			case flag_10bit::INFLOW: // assuming only x-direction
				U.at(i).at(0)=1;
				V.at(i).at(0) = 0;
				break;

			default:
				break;
		}
	}

	// ---------------------------------------------- Top-Left Corner ----------------------------------------------------------------------
	switch(flag.at(0).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
			if (B_N(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_N at top-left boundary "));
			}

			if (B_W(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_W at top-left boundary "));
			}

			if (B_S(flag.at(0).at(jmax+1)) == 1){
				V.at(0).at(jmax) = 0;
				U.at(0).at(jmax+1) = -U.at(0).at(jmax);
			}

			if (B_E(flag.at(0).at(jmax+1)) == 1){
				U.at(0).at(jmax+1) = 0;
				V.at(0).at(jmax) = -V.at(1).at(jmax);
				V.at(0).at(jmax+1) = -V.at(1).at(jmax+1);
			}

			if (B_NE(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary "));
			}

			if (B_NW(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary "));
			}

			if (B_SE(flag.at(0).at(jmax+1)) == 1){
				U.at(0).at(jmax+1)=0;
				V.at(0).at(jmax)=0;
				V.at(0).at(jmax+1) = -V.at(1).at(jmax+1);
			}

			if (B_SW(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary "));
			}
			break;

		case flag_10bit::FREE_SLIP://Free Slip conditions
			if (B_E(flag.at(0).at(jmax+1)) == 1){
				U.at(0).at(jmax+1) = 0;
				V.at(0).at(jmax+1) = V.at(1).at(jmax+1);
				V.at(0).at(jmax) = V.at(1).at(jmax);
			}

			if (B_W(flag.at(0).at(jmax+1)) == 1) {
				throw std::runtime_error(std::string("impossible to have B_W at top-left boundary "));
			}

			if (B_N(flag.at(0).at(jmax+1)) == 1) {
				throw std::runtime_error(std::string("impossible to have B_N at top-left boundary "));
			}

			if (B_S(flag.at(0).at(jmax+1)) == 1){
				V.at(0).at(jmax) = 0;
				U.at(0).at(jmax+1) = U.at(0).at(jmax);
			}

			if (B_NE(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary "));
			}

			if (B_NW(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary "));
			}

			if (B_SE(flag.at(0).at(jmax+1)) == 1){
				U.at(0).at(jmax+1)=0;
				V.at(0).at(jmax)=0;
				V.at(0).at(jmax+1) = V.at(1).at(jmax+1);
			}

			if (B_SW(flag.at(0).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary "));
			}
			break;

		case flag_10bit::OUTFLOW: // assuming only x-direction
			throw std::runtime_error(std::string("impossible to have outflow at top-left boundary (only x direction assumed"));
			break;

		case flag_10bit::INFLOW: // assuming only x-direction
			U.at(0).at(jmax+1)=1;
			V.at(0).at(jmax+1) = 0;
			V.at(0).at(jmax) = 0;
			break;

		default:
			break;
	}

	// ------------------------------------------------- Bottom-left Corner ------------------------------------------------------------------
	switch(flag.at(0).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
			if (B_N(flag.at(0).at(0)) == 1){
				V.at(0).at(0) = 0;
				U.at(0).at(0) = -U.at(0).at(1);
			}

			if (B_W(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary"));
			}

			if (B_S(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary"));
			}

			if (B_E(flag.at(0).at(0)) == 1){
				U.at(0).at(0) = 0;
				V.at(0).at(0) = -V.at(1).at(0);
			}

			if (B_NE(flag.at(0).at(0)) == 1){
				U.at(0).at(0) = 0;
				V.at(0).at(0) = 0;
			}

			if (B_NW(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary"));
			}

			if (B_SE(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary"));
			}

			if (B_SW(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary"));
			}
			break;

		case flag_10bit::FREE_SLIP://Free Slip conditions
			if (B_E(flag.at(0).at(0)) == 1){
				U.at(0).at(0) = 0;
				V.at(0).at(0) = V.at(1).at(0);
			}

			if (B_W(flag.at(0).at(0)) == 1) {
				throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary"));
			}

			if (B_N(flag.at(0).at(0)) == 1) {
				V.at(0).at(0) = 0;
				U.at(0).at(0) = U.at(0).at(1);
			}

			if (B_S(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary"));
			}

			if (B_NE(flag.at(0).at(0)) == 1){
				U.at(0).at(0) = 0;
				V.at(0).at(0) = 0;
			}

			if (B_NW(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary"));
			}

			if (B_SE(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary"));
			}

			if (B_SW(flag.at(0).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary"));
			}
			break;

		case flag_10bit::OUTFLOW: // assuming only x-direction
			throw std::runtime_error(std::string("impossible to have outflow at bottom-left boundary (only x direction)"));
			break;

		case flag_10bit::INFLOW: // assuming only x-direction
			U.at(0).at(0)= 1;
			V.at(0).at(0) = 0;
			break;

		default:
			break;
	}

	// ----------------------------------------------------- Top-right corner -----------------------------------------------------------
	switch(flag.at(imax+1).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
			if (B_N(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_N at top-right boundary"));
			}

			if (B_W(flag.at(imax+1).at(jmax+1)) == 1){
				U.at(imax).at(jmax+1) = 0;
				V.at(imax+1).at(jmax) = -V.at(imax).at(jmax);
				V.at(imax+1).at(jmax+1) = -V.at(imax).at(jmax+1);
			}

			if (B_S(flag.at(imax+1).at(jmax+1)) == 1){
				V.at(imax+1).at(jmax) = 0;
				U.at(imax).at(jmax+1) = -U.at(imax).at(jmax);
				U.at(imax+1).at(jmax+1) = -U.at(imax+1).at(jmax);
			}

			if (B_E(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_E at top-right boundary"));
			}

			if (B_NE(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at top-right boundary"));
			}

			if (B_NW(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary"));
			}

			if (B_SE(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary"));
			}

			if (B_SW(flag.at(imax+1).at(jmax+1)) == 1){
				U.at(imax).at(jmax+1) = 0;
				U.at(imax+1).at(jmax+1) = -U.at(imax+1).at(jmax);
				V.at(imax+1).at(jmax) = 0;
				V.at(imax+1).at(jmax+1) = -V.at(imax).at(jmax+1);
			}
			break;

		case flag_10bit::FREE_SLIP://Free Slip conditions
			if (B_E(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_E at top-right boundary"));
			}

			if (B_W(flag.at(imax+1).at(jmax+1)) == 1) {
				U.at(imax).at(jmax+1) = 0;
				V.at(imax+1).at(jmax) = V.at(imax).at(jmax);
				V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
			}

			if (B_N(flag.at(imax+1).at(jmax+1)) == 1) {
				throw std::runtime_error(std::string("impossible to have B_N at top-right boundary"));
			}

			if (B_S(flag.at(imax+1).at(jmax+1)) == 1){
				V.at(imax+1).at(jmax) = 0;
				U.at(imax).at(jmax+1) = U.at(imax).at(jmax);
				U.at(imax+1).at(jmax+1) = U.at(imax+1).at(jmax);
			}

			if (B_NE(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at top-right boundary"));
			}

			if (B_NW(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary"));
			}

			if (B_SE(flag.at(imax+1).at(jmax+1)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary"));
			}

			if (B_SW(flag.at(imax+1).at(jmax+1)) == 1){
				U.at(imax).at(jmax+1) = 0;
				U.at(imax+1).at(jmax+1) = U.at(imax+1).at(jmax);
				V.at(imax+1).at(jmax) = 0;
				V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
			}
			break;

		case flag_10bit::OUTFLOW: // assuming only x-direction
			U.at(imax+1).at(jmax+1) = U.at(imax).at(jmax+1);
			V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
			V.at(imax+1).at(jmax) = V.at(imax).at(jmax);
			break;

		case flag_10bit::INFLOW: // assuming only x-direction
			U.at(imax+1).at(jmax+1)=1;
			V.at(imax+1).at(jmax+1) = 0;
			V.at(imax+1).at(jmax) = 0;
			break;

		default:
			break;
	}

	// ------------------------------------------------------ Bottom-right corner -----------------------------------------------------------
	switch(flag.at(imax+1).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
			if (B_N(flag.at(imax+1).at(0)) == 1){
				V.at(imax+1).at(0) = 0;
				U.at(imax).at(0) = -U.at(imax).at(1);
				U.at(imax+1).at(0) = -U.at(imax+1).at(1);
			}

			if (B_W(flag.at(imax+1).at(0)) == 1){
				U.at(imax).at(0) = 0;
				V.at(imax+1).at(0) = -V.at(imax).at(0);
			}

			if (B_S(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary"));
			}

			if (B_E(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary"));
			}

			if (B_NE(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary"));
			}

			if (B_NW(flag.at(imax+1).at(0)) == 1){
				U.at(imax).at(0) = 0;
				U.at(imax+1).at(0) = - U.at(imax+1).at(1);
				V.at(imax+1).at(0) = 0;
			}

			if (B_SE(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary"));
			}

			if (B_SW(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary"));
			}
			break;

		case flag_10bit::FREE_SLIP://Free Slip conditions
			if (B_E(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary"));
			}

			if (B_W(flag.at(imax+1).at(0)) == 1) {
				U.at(imax).at(0) = 0;
				V.at(imax+1).at(0) = V.at(imax).at(0);
			}

			if (B_N(flag.at(imax+1).at(0)) == 1) {
				V.at(imax+1).at(0) = 0;
				U.at(imax).at(0) = U.at(imax).at(1);
				U.at(imax+1).at(0) = U.at(imax+1).at(1);
			}

			if (B_S(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary"));
			}

			if (B_NE(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary"));
			}

			if (B_NW(flag.at(imax+1).at(0)) == 1){
				U.at(imax).at(0) = 0;
				U.at(imax+1).at(0) = U.at(imax+1).at(1);
				V.at(imax+1).at(0) = 0;
			}

			if (B_SE(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary"));
			}

			if (B_SW(flag.at(imax+1).at(0)) == 1){
				throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary"));
			}
			break;

		case flag_10bit::OUTFLOW: // assuming only x-direction
			U.at(imax+1).at(0) = U.at(imax).at(0);
			V.at(imax+1).at(0) = V.at(imax).at(0);
			break;

		case flag_10bit::INFLOW: // assuming only x-direction
			U.at(imax+1).at(0)=1;
			V.at(imax+1).at(0) = 0;
			break;

		default:
			break;
	}

	// ---------------------------------------------------------- Inner Cells --------------------------------------------------------------------
	for(int i = 1; i<=imax; ++i){

		for(int j = 1; j<=jmax; ++j){

			switch(flag.at(i).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW )){
				case flag_10bit::NO_SLIP: //No slip conditions
					if (B_N(flag.at(i).at(j)) == 1){
						V.at(i).at(j) = 0;
						U.at(i-1).at(j) = -U.at(i-1).at(j+1);
						U.at(i).at(j) = -U.at(i).at(j+1);
					}

					if (B_W(flag.at(i).at(j)) == 1){
						U.at(i-1).at(j) = 0;
						V.at(i).at(j-1) = -V.at(i-1).at(j-1);
						V.at(i).at(j) = -V.at(i-1).at(j);
					}

					if (B_S(flag.at(i).at(j)) == 1){
						V.at(i).at(j-1) = 0;
						U.at(i-1).at(j) = -U.at(i-1).at(j-1);
						U.at(i).at(j) = -U.at(i).at(j-1);
					}

					if (B_E(flag.at(i).at(j)) == 1){
						U.at(i).at(j) = 0;
						V.at(i).at(j-1) = -V.at(i+1).at(j-1);
						V.at(i).at(j) = -V.at(i+1).at(j);
					}

					if (B_NE(flag.at(i).at(j)) == 1){
						U.at(i).at(j) = 0;
						U.at(i-1).at(j) = -U.at(i-1).at(j+1);
						V.at(i).at(j) = 0;
						V.at(i).at(j-1) = -V.at(i+1).at(j-1);
					}

					if (B_NW(flag.at(i).at(j)) == 1){
						U.at(i-1).at(j) = 0;
						U.at(i).at(j) = - U.at(i).at(j+1);
						V.at(i).at(j) = 0;
						V.at(i).at(j-1) = -V.at(i-1).at(j-1);
					}

					if (B_SE(flag.at(i).at(j)) == 1){
						U.at(i).at(j)=0;
						U.at(i-1).at(j) = -U.at(i-1).at(j-1);
						V.at(i).at(j-1)=0;
						V.at(i).at(j) = -V.at(i+1).at(j);
					}

					if (B_SW(flag.at(i).at(j)) == 1){
						U.at(i-1).at(j) = 0;
						U.at(i).at(j) = -U.at(i).at(j-1);
						V.at(i).at(j-1) = 0;
						V.at(i).at(j) = -V.at(i-1).at(j);
					}
					break;

				case flag_10bit::FREE_SLIP://Free Slip conditions
					if (B_E(flag.at(i).at(j)) == 1){
						U.at(i).at(j) = 0;
						V.at(i).at(j) = V.at(i+1).at(j);
						V.at(i).at(j-1) = V.at(i+1).at(j-1);
					}

					if (B_W(flag.at(i).at(j)) == 1) {
						U.at(i-1).at(j) = 0;
						V.at(i).at(j-1) = V.at(i-1).at(j-1);
						V.at(i).at(j) = V.at(i-1).at(j);
					}

					if (B_N(flag.at(i).at(j)) == 1) {
						V.at(i).at(j) = 0;
						U.at(i-1).at(j) = U.at(i-1).at(j+1);
						U.at(i).at(j) = U.at(i).at(j+1);
					}

					if (B_S(flag.at(i).at(j)) == 1){
						V.at(i).at(j-1) = 0;
						U.at(i-1).at(j) = U.at(i-1).at(j-1);
						U.at(i).at(j) = U.at(i).at(j-1);
					}

					if (B_NE(flag.at(i).at(j)) == 1){
						U.at(i).at(j) = 0;
						U.at(i-1).at(j) = U.at(i-1).at(j+1);
						V.at(i).at(j) = 0;
						V.at(i).at(j-1) = V.at(i+1).at(j-1);
					}

					if (B_NW(flag.at(i).at(j)) == 1){
						U.at(i-1).at(j) = 0;
						U.at(i).at(j) = U.at(i).at(j+1);
						V.at(i).at(j) = 0;
						V.at(i).at(j-1) = V.at(i-1).at(j-1);
					}

					if (B_SE(flag.at(i).at(j)) == 1){
						U.at(i).at(j)=0;
						U.at(i-1).at(j) = U.at(i-1).at(j-1);
						V.at(i).at(j-1)=0;
						V.at(i).at(j) = V.at(i+1).at(j);
					}

					if (B_SW(flag.at(i).at(j)) == 1){
						U.at(i-1).at(j) = 0;
						U.at(i).at(j) = U.at(i).at(j-1);
						V.at(i).at(j-1) = 0;
						V.at(i).at(j) = V.at(i-1).at(j);
					}
					break;

				case flag_10bit::OUTFLOW: // assuming only x-direction
					U.at(i).at(j) = U.at(i-1).at(j);
					V.at(i).at(j) = V.at(i-1).at(j);
					V.at(i).at(j-1) = V.at(i-1).at(j-1);
					break;

				case flag_10bit::INFLOW: // assuming only x-direction
					U.at(i).at(j)=1;
					V.at(i).at(j) = 0;
					V.at(i).at(j-1) = 0;
					break;

				default:
					break;
			}
		}
	}

}

// ------------------------------------------- check conditions for obstacle cells -----------------------------------------------------

int B_E(unsigned int flag){
    return (((flag & flag_10bit::E) == flag_10bit::E) && !(((flag & flag_10bit::N) == flag_10bit::N) || ((flag & flag_10bit::S) == flag_10bit::S)));
}

int B_W(unsigned int flag){
    return (((flag & flag_10bit::W) == flag_10bit::W) && !(((flag & flag_10bit::N) == flag_10bit::N) || ((flag & flag_10bit::S) == flag_10bit::S)));
}

int B_N(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) && !(((flag & flag_10bit::E) == flag_10bit::E) || ((flag & flag_10bit::W) == flag_10bit::W)));
}

int B_S(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) && !(((flag & flag_10bit::E) == flag_10bit::E) || ((flag & flag_10bit::S) == flag_10bit::W)));
}

int B_NE(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) &&  ((flag & flag_10bit::E) == flag_10bit::E));
}

int B_NW(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) &&  ((flag & flag_10bit::W) == flag_10bit::W));
}

int B_SE(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) &&  ((flag & flag_10bit::E) == flag_10bit::E));
}

int B_SW(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) &&  ((flag & flag_10bit::W) == flag_10bit::W));
}
