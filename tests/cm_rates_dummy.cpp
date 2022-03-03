#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	//TestModel M;
	Environment E;
	Species<TestModel> spp;

	Solver S(SOLVER_CM);
	S.use_log_densities = false;  // FIXME:: Make sure that oder of calling this before/after reset state doesnt matter
	S.control.cm_grad_dx = 0.001;
	S.addSpecies(25, 0, 1, false, &spp, 4);
	S.resetState();
	S.initialize();
	S.print();
	
	E.computeEnv(0,&S, S.state.begin(), S.rates.begin());
	//cout << E.evalEnv(0,0) << endl;

	S.setEnvironment(&E);
	S.calcRates_CM(1, S.state.begin(), S.rates.begin());  // dummy rates calc rates(X=X0, U=U0, t=1, E=E(U0))
	//S.print();
	//	S.step_to(1);

	//vector <double> rates_exp = {
		  //0.24871666,  0.24831871,  0.24712487,  0.24513514,  0.24234951,  0.23876799,  0.23439058,  0.22921727,
		  //0.22324807,  0.21648298,  0.20892199,  0.20056511,  0.19141234,  0.18146367,  0.17071911,  0.15917866,
		  //0.14684231,  0.13371007,  0.11978194,  0.10505792,  0.08953800,  0.07322218,  0.05611048,  0.03820288,
		  //0.01949939,  0.00000000, -3.07367787, -2.48964104, -2.02468490, -1.65220631, -1.35210765, -1.10906811,
		 //-0.91130994, -0.74970832, -0.61714231, -0.50801705, -0.41790887, -0.34329937, -0.28137461, -0.22987235,
		 //-0.18696499, -0.15116931, -0.12127634, -0.09629670, -0.07541760, -0.05796891, -0.04339629, -0.03123970,
		 //-0.02111622, -0.01270624, -0.00574230, 0.00000000 };

	vector <double> rates_exp = {
		0,0,
		0.01949939,-0.0057423,
		0.03820288,-0.01270624,
		0.05611048,-0.02111622,
		0.07322218,-0.0312397,
		0.089538,-0.04339629,
		0.10505792,-0.05796891,
		0.11978194,-0.0754176,
		0.13371007,-0.0962967,
		0.14684231,-0.12127634,
		0.15917866,-0.15116931,
		0.17071911,-0.18696499,
		0.18146367,-0.22987235,
		0.19141234,-0.28137461,
		0.20056511,-0.34329937,
		0.20892199,-0.41790887,
		0.21648298,-0.50801705,
		0.22324807,-0.61714231,
		0.22921727,-0.74970832,
		0.23439058,-0.91130994,
		0.23876799,-1.10906811,
		0.24234951,-1.35210765,
		0.24513514,-1.65220631,
		0.24712487,-2.0246849,
		0.24831871,-2.48964104,
		0.24871666,-3.07367787,
		-2,0,-20,-30,
		-2,0,-19.2,-28.8,
		-2,0,-18.4,-27.6,
		-2,0,-17.6,-26.4,
		-2,0,-16.8,-25.2,
		-2,0,-16,-24,
		-2,0,-15.2,-22.8,
		-2,0,-14.4,-21.6,
		-2,0,-13.6,-20.4,
		-2,0,-12.8,-19.2,
		-2,0,-12,-18,
		-2,0,-11.2,-16.8,
		-2,0,-10.4,-15.6,
		-2,0,-9.6,-14.4,
		-2,0,-8.8,-13.2,
		-2,0,-8,-12,
		-2,0,-7.2,-10.8,
		-2,0,-6.4,-9.6,
		-2,0,-5.6,-8.4,
		-2,0,-4.8,-7.2,
		-2,0,-4,-6,
		-2,0,-3.2,-4.8,
		-2,0,-2.4,-3.6,
		-2,0,-1.6,-2.4,
		-2,0,-0.8,-1.2,
		-2,0,0,0};

	
	for (int i=0; i< S.rates.size(); ++i){
		cout << S.rates[i] << " " << rates_exp[i] << endl;
		if ( fabs(S.rates[i] - rates_exp[i]) > 1e-6) return 1;
	}
	
	return 0;
	
}


