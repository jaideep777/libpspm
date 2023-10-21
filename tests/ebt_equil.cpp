#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

std::vector <double> diff(vector <double> breaks){
	std::vector<double> mids(breaks.size()-1);
	for (size_t i=0; i<mids.size(); ++i) mids[i] = (breaks[i]+breaks[i+1])/2;
	return mids;
}

int main(){

	// TestModel M;
	// Species<TestModel> spp(M);
	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_EBT);
	S.setEnvironment(&E);
	S.addSpecies({25}, {0}, {1}, {false}, &spp, 4, -1);

	S.print();
	
	E.computeEnv(0, &S, S.state.begin(), S.rates.begin());
	cout << E.evalEnv(0,0) << endl;
	if (fabs(E.evalEnv(0,0) - 0.3802) > 1e-5) return 1;

	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	
	ofstream fout("ebt_testmodel_equil.txt");

	vector<double> breaks = myseq(0,1,26);
	vector<double> mids = diff(breaks);

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";

		// vector<double> breaks = myseq(0,1,26);
		// vector<double> v = S.getDensitySpecies(0, breaks, Spline::QUADRATIC);
		// for (auto y : v) fout << y << "\t";
		// fout << endl;
	}

	fout.close();

	S.print();	
	cout << setprecision(6) << S.u0_out(S.current_time)[0] << endl;
	if (abs(S.u0_out(S.current_time)[0]-0.97504) < 2e-5) return 0;
	else return 1;
}

