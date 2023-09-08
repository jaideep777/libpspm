#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"

using namespace std;

#include "red.h"



inline std::vector <double> logseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = exp(log(from) + i*(log(to)-log(from))/(len-1));
	return x;
}


int main(){

	Species<RED_Plant> spp;
	LightEnvironment E;

	Solver S(SOLVER_EBTN);
	S.control.ebt_ucut = 1e-20;

	S.setEnvironment(&E);
	std::vector<double> xb;
	xb.push_back(1);
	xb.push_back(0.5);
	std::vector<double> xm;
	xm.push_back(1e4);
	xm.push_back(100);
	std::vector<bool> logBreaks;
	logBreaks.push_back(true);
	logBreaks.push_back(false);

	S.addSpecies(10, xb, xm, logBreaks, &spp, 0);
	
	
	S.resetState();
	S.initialize();

	S.printODEmethod();
	// S.print();
	
	
	ofstream fout("ebtn_Redmodel.txt");


	std::ofstream cohortprint;
	cohortprint.open(std::string("cohort_vector_ebtn_Redmodel.csv").c_str());

	for (double t=0; t <= 10; t=t+1) {
		S.step_to(t);
		S.printCohortVector(cohortprint);
		cout << "Finished step to function for t = " << t << std::endl;
		// S.printCohortVector(std::cout);

		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
		cout << S.current_time << " " << S.species_vec[0]->xsize() << "\n";
		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
		
		// vector <double> dist = S.getDensitySpecies(0, logseq(1, 1e6, 150));
		// for (auto y : dist) fout << y << "\t";
		fout << endl;
	}

	fout.close();
	cohortprint.close();

	// // expected 38.1128953 (numerical R), 37.5845 (analytical)
	// cout << S.newborns_out(5000)[0] << endl; 
	// //if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	// //else return 1;

}

