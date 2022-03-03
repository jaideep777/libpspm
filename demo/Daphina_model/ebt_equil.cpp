#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "daphnia.h"

int main(){

	Species<Daphnia> spp;
	Environment E;

	Solver S(SOLVER_EBT);

	S.addSpecies(100, 0, 1, false, &spp, 0, -1);
	S.addSystemVariables(1);  // this can be done either before or after addSpecies()

	S.control.ebt_ucut = 1e-20;

	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	S.state[0] = E.K;
	//S.print();
	
	
	ofstream fout("ebt_Daphnia.txt");

	for (double t=0.05; t <= 100; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out(t)[0] << "\t";
		cout << S.current_time << " " << S.state[0] << "\n";
		
		vector <double> dist = S.getDensitySpecies_EBT(0, 30);
		for (auto y : dist) fout << y << "\t";
		
		fout << endl;
	}

	fout.close();

	// expected 38.1128953 (numerical R), 37.5845 (analytical)
	cout << S.newborns_out(100)[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

