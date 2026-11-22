#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
#include "individual_base.h"
using namespace std;

#define _FISHERY_DT_ (1/12.0)

#include "cod_model.h"

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

int main(){

	Species<PSPMFish> spp;

	SeaEnv E;

	Solver S(SOLVER_IEBT);
	S.control.ode_ifmu_stepsize = _FISHERY_DT_;
	S.control.cohort_insertion_dt = 1;
	S.setEnvironment(&E);
	S.addSpecies({30}, {0.5}, {30.5}, {false}, &spp, 0, -1);
	S.species_vec[0]->set_xb({0.001});
	S.print();

	ofstream fout1;
	fout1.open("iebt_fish_u.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	// for (int i=0; i<50; ++i){
	// 	S.copyCohortsToState();
	// 	S.stepU_iFMU(0, S.state, S.rates, 0.02);
	// 	S.copyStateToCohorts(S.state.begin());
	// 	// S.print();
	// }

	ofstream foutx("iebt_fish_x.txt");

	ofstream fout("iebt_fishmodel_equil.txt");
	ofstream foutp("iebt_fishmodel_communityprops.txt");

	for (double t=0; t <= 100; t=t+1) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << " " << static_cast<SeaEnv*>(S.env)->ssb/1e9 << " " << S.species_vec[0]->xsize() << "\n";

		// S.print();
		// // cout.flush();
		// cout << "U integral = " << S.state_integral([](int i, double t){return 1;}, t, 0) << "\n";

		// vector<double> breaks = myseq(0,10,201);
		// vector<double> vx = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		// for (auto x : vx) foutx << x << "\t";
		// foutx << endl;

		for (double a=1; a<30; a+=1){
			double U_mat = S.integrate_wudx_above([&S](int i, double t){return (S.species_vec[0]->getC(i) == 1)? 1:0;}, 0, {a}, 0)
						- S.integrate_wudx_above([&S](int i, double t){return (S.species_vec[0]->getC(i) == 1)? 1:0;}, 0, {a-0.99999}, 0);
			double U_tot = S.integrate_wudx_above([&S](int i, double t){return 1;}, 0, {a}, 0)
						- S.integrate_wudx_above([&S](int i, double t){return 1;}, 0, {a-0.99999}, 0);
			double mat = (U_mat)/(U_tot+1e-6);
			// cout << "Umat/Utot = " << mat << " at x = " << a << "\n";
			foutp << t << "\t" << a-1 << "\t" << "U_tot\t" << U_tot << "\n";
			foutp << t << "\t" << a-1 << "\t" << "mat\t" << mat << "\n";
		}
		foutp << t << "\t" << "NA" << "\t" << "ssb\t" << static_cast<SeaEnv*>(S.env)->ssb << "\n";
		foutp << t << "\t" << "NA" << "\t" << "tsb\t" << static_cast<SeaEnv*>(S.env)->tsb << "\n";
	}

	fout1.open("iebt_fish_u1.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	fout.close();
	foutx.close();
	foutp.close();

	// S.print();

	// cout << S.u0_out(S.current_time)[0] << endl; 
	// if (abs(S.u0_out(S.current_time)[0] - 1.30639) < 1e-5) return 0;
	// else return 1;

}

