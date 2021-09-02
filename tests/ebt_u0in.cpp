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

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_EBT);
	S.addSpecies(25, 0, 1, false, &spp, 4, 2);
	S.species_vec[0]->set_bfin_is_u0in(true);
	S.resetState();
	S.initialize();
	S.print();
	
	E.computeEnv(0, &S);
	cout << E.evalEnv(0,0) << endl;

	S.setEnvironment(&E);
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	
	ofstream fout("ebt_testmodel.txt");

	vector<double> breaks = myseq(0,1,26);
	vector<double> mids = diff(breaks);

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out()[0] << "\t";

		//cout << "HERE" << endl;	
		vector<double> v = S.getDensitySpecies_EBT(0, 26);
				
		cout << S.current_time << " " << S.species_vec[0]->xsize() << " " << S.u0_out()[0] << endl;
		//for (auto y : S.state) fout << y << "\t";
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	S.print();	
	cout << S.u0_out()[0] << endl;
	if (abs(S.u0_out()[0]-1.436407) < 2e-5) return 0;
	else return 1;
}

