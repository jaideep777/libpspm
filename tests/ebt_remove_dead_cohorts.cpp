#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

template<class T>
bool almostEqual(const vector<T> &a, const vector <T> &b){
	if (a.size() != b.size()) return false;
	return std::equal(a.begin(), a.end(), b.begin(), [](const T& x, const T& y){return abs(x-y)<1e-5;});
}


int main(){
	TestModel M;
	Environment E;

	Solver<TestModel,Environment> S(SOLVER_EBT);
	S.use_log_densities = false;
	//S.control.cm_grad_dx = 0.001;
	S.addSpecies(5, 0, 1, false, &M, {"mort", "vs", "ha", "sa"}, 2);
	S.addSpecies(8, 0, 1, false, &M, {"mort", "vs", "ha", "sa"}, 2);
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;
	
	auto itu1 = S.get_species(0)->get_iterators(S.state).get("u");	
	auto itu2 = S.get_species(1)->get_iterators(S.state).get("u");

	*(itu1+2) = *(itu1+3) = 1e-11;
	*(itu2+4) = *(itu2+6) = *(itu2+7) = 1e-11;

	S.print();

	S.removeDeadCohorts_EBT();

	S.print();

	return 0;
	
}

