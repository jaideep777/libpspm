#include <cassert>
#include <fstream>
#include "solver.h"
#include "index_utils.h"
#include "linear_system.h"
using namespace std;


// Note: This IFMU solver works only for non-negative growth functions
void Solver::stepU_iFMU(double t, vector<double> &S, vector<double> &dSdt, double dt){
	
	// static matrixes for solving MU=Z
	static Matrix M;
	static vector<double> Z;

	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *U = &(*its); // Since FMU only has U in state, start of species is actually U
		int J = spp->J;			// xsize of species.
		vector<vector<double>> &x = spp->x;
		vector<vector<double>> &X = spp->X;
		vector<vector<double>> &h = spp->h;
		
		vector <vector<double>> growthArray(J);
		for (int i=0; i<J; ++i) growthArray[i] = spp->growthRate(i, t, env);
		
		double birthFlux;
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			birthFlux = spp->birth_flux_in * pe;
		}
		//cout << t << "\t | B*pe = " << calcSpeciesBirthFlux(s,t) << " * " << pe << " -> " << birthFlux << "\n";

		// FIXME: Eventually replace this with a standard sparse matrix library
		// Resize matrix to correct dims (J x J) and fill it with zeros
		M.clear();
		M.resize(J);
		for (auto &m : M) m.resize(J, 0);
		Z.clear(); Z.resize(J, 0);

		// handle boundary condition at corner cell
		double dV = 1; 
		for (int ih=0; ih<h.size(); ++ih) dV *= h[ih][0]; // calc volume of corner cell
		double A0 = 1 + spp->mortalityRate(0,t,env)*dt;
		for (int m=0; m<spp->istate_size; ++m) A0 += growthArray[0][m]*dt/h[m][0];
		M[0][0] = A0;
		Z[0] = U[0] + birthFlux/dV*dt;
		// cout << "B = " << birthFlux << '\n';

		// traverse remaining cells to fill matrix
		for (int j=1; j<J; ++j){ // loop over all cohorts except the corner cohort
			vector<int> id = id_utils::index(j, spp->dim_centres); // get cell index. Remember that cells are labelled by upper corner
			vector<double> Xcell = id_utils::coord_value(id, X); // upper corner value of cell at id

			// Calculate diagonal elements Aj 
			double Aj = 1 + spp->mortalityRate(j,t,env)*dt;
			for (int m=0; m<spp->istate_size; ++m){
				// use backward difference if growth rate is +ve OR if cell is on the upper edge wrt axis m
				if (growthArray[j][m] >= 0 || 
					id[m] == spp->dim_centres[m]-1 ){ 
					Aj += growthArray[j][m]*dt/h[m][id[m]];
				}
				else{ // else use forward difference when growth rate is negative
					Aj -= growthArray[j][m]*dt/h[m][id[m]+1];
				}
			}
			M[j][j] = Aj;
			Z[j] = U[j];

			// Calculate off-diagonal elements
			for (int m=0; m<spp->istate_size; ++m){
				// use backward difference if growth rate is +ve OR if cell is on the upper edge wrt axis m
				if (growthArray[j][m] >= 0 || 
					id[m] == spp->dim_centres[m]-1 ){ 
					if (id[m] == 0) continue; // This term is 0 if id[m] = 0, so do nothing.
					
					int j_minus_1m = id_utils::loc_minus1k(id, m, spp->dim_centres);
					M[j][j_minus_1m] = -growthArray[j_minus_1m][m]*dt/h[m][id[m]];
				}
				else{ // else use forward difference when growth rate is negative
					int j_plus_1m = id_utils::loc_plus1k(id, m, spp->dim_centres);
					M[j][j_plus_1m] = growthArray[j_plus_1m][m]*dt/h[m][id[m]+1];
				}
			}
		}

		// // ~~~ Old method ~~~
		// // Vector Unew = lupSolve(M, Z); // FIXME: Replace with Eigen or suchlike
		// // for (int i=0; i<J; ++i) U[i] = Unew[i];
		// double B0  = 1 + dt/h[0][0]*growthArray[0][0] + dt*spp->mortalityRate(0, t, env);
		// double C0 = spp->getU(0) + dt/h[0][0]*birthFlux;
		// U[0] = C0/B0;

		// // O1 scheme
		// for (int w = 1; w < J; ++w){
		// 	double Aw = -growthArray[w-1][0]*dt/h[0][w];
		// 	double Bw  = 1 + dt/h[0][w]*growthArray[w][0] + dt*spp->mortalityRate(w, t, env);
		// 	double Cw = spp->getU(w);

		// 	U[w] = (Cw - Aw*U[w-1])/Bw;
		// }

		// // *** O2 scheme ****
		// // O1 scheme for w = 0,1
		// for (int w = 1; w < 2; ++w){
		// 	double Aw = -growthArray[w-1]*dt/h[w];
		// 	double Bw  = 1 + dt/h[w]*growthArray[w] + dt*spp->mortalityRate(w, spp->getX(w), t, env);
		// 	double Cw = spp->getU(w);

		// 	U[w] = (Cw - Aw*U[w-1])/Bw;
		// }
		// for (int w = 2; w < J; ++w){
				
		// 	double Mw = spp->mortalityRate(w, spp->getX(w), t, env);
			
		// 	// O1 scheme
		// 	double Aw = -growthArray[w-1]*dt/h[w];
		// 	double Bw  = 1 + dt/h[w]*growthArray[w] + dt*Mw;
		// 	double Cw = spp->getU(w);

		// 	double Uw1 = (Cw - Aw*U[w-1])/Bw;

		// 	// O2 scheme
		// 	double A2w = 1 + 3*dt/(2*h[w])*growthArray[w] + dt*Mw;
		// 	double B2w = -4*growthArray[w-1]*dt/(2*h[w]);
		// 	double C2w = growthArray[w-2]*dt/(2*h[w]);
		// 	double D2w = spp->getU(w);

		// 	double Uw2 = (D2w - B2w*U[w-1] - C2w*U[w-2])/A2w;

		// 	double phi = control.ifmu_order - 1; // 0 for O1, 1 for O2, or in-between
 		// 	//if (Uw2 < 0) phi = 0; // use O1 if O2 density is negative. FIXME: this condition causes wrong results in Daphnia model!   
		// 	U[w] = phi*Uw2 + (1-phi)*Uw1;
		// }


		its += J*(1+spp->n_accumulators);

	}
	
}


