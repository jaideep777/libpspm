#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
#include "individual_base.h"
using namespace std;


class WaveEnv : public EnvironmentBase{
	public:
	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
	}

};

class Wave : public IndividualBase<2>{
	public:
	double init_density(void * _env, double bf){
		return exp(-10*(x[0]-2)*(x[0]-2) - 10*(x[1]-4)*(x[1]-4));
		// return (x[0] > 2 && x[0] < 3 && x[1] > 3 && x[1] < 4)? 1 : 0;
	}

	void set_size(const array<double,2>& x){
	}

	array<double,2> growthRate(double t, void * env){
		return {1,2};
	}
	double mortalityRate(double t, void * env){
		return 0;
	}
	double birthRate(double t, void * env){
		return 0.1;
	}

};



std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

int main(){

	Species<Wave> spp;
	WaveEnv E;

	Solver S(SOLVER_IFMU);
	S.control.ode_ifmu_stepsize = 0.02;
	S.setEnvironment(&E);
	S.addSpecies({25, 25}, {0, 0}, {10,10}, {false, false}, &spp, 0, -1);
	S.print();
	ofstream fout1;
	fout1.open("ifmu2d_u.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	// for (int i=0; i<50; ++i){
	// 	S.copyCohortsToState();
	// 	S.stepU_iFMU(0, S.state, S.rates, 0.02);
	// 	S.copyStateToCohorts(S.state.begin());
	// 	// S.print();
	// }

	ofstream fout("ifmu_testmodel_equil.txt");

	for (double t=0; t <= 1; t=t+1) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
		
		// vector<double> breaks = myseq(0,1,26);
		// vector<double> v = S.getDensitySpecies1D(0, 0, breaks, Spline::QUADRATIC);
		// for (auto y : v) fout << y << "\t";
		// fout << endl;
	}
	fout1.open("ifmu2d_u1.txt");
	S.species_vec[0]->printCohortVector(fout1);
	fout1.close();

	fout.close();

	// cout << S.u0_out(S.current_time)[0] << endl; 
	// if (abs(S.u0_out(S.current_time)[0] - 1.30639) < 1e-5) return 0;
	// else return 1;

}

