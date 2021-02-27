#include <solver.h>
#include "test_model_2_ms.h"

int main(){
	TestModel M;
	Environment E;

	Solver<TestModel,Environment> S(SOLVER_CM);
	S.use_log_densities = true;
	
	//Solver<TestModel> S({0,1,2,3}, SOLVER_CM);
	S.addSpecies({0,1,2,3}, &M, {}, 2);
	if(S.get_species(0)->xb != 0) return 1;

	S.get_species(0)->set_bfin_is_u0in(true);
	S.get_species(0)->set_inputBirthFlux(2);
	S.resetState();
	//S.initialize();
	//S.print();

	S.state = {0,1,2,3, log(1),log(2),log(4),log(8)};
	S.get_species(0)->u0_save = S.get_u0(0,0);
	//S.print();

	auto w = [](double x, double t){return 1;};

	cout << S.integrate_x(w, 0, S.state, 0) << " " << S.integrate_wudx_above(w, 0, 0, S.state, 0) << endl;
	
	if(abs(S.integrate_x(w, 0, S.state, 0) - 10.5) > 1e-5) return 1;
	if(abs(S.integrate_x(w, 0, S.state, 0) - S.integrate_wudx_above(w, 0, 0, S.state, 0)) > 1e-5) return 1;
	//if(abs(S.integrate_wudx_above(w, 0, 1.5, S.state) - 7.75) < 2e-5) return 1;


	S.state = {1,2,3,4, log(1),log(2),log(4),log(8)};
	cout << S.integrate_x(w, 0, S.state, 0) << " " << S.integrate_wudx_above(w, 0, 1.1, S.state, 0) << endl;
	if(abs(S.integrate_x(w, 0, S.state, 0) - 12) > 1e-5) return 1;
	if(abs(S.integrate_wudx_above(w, 0, 1.1, S.state, 0) - 10.5) > 1e-5) return 1;


	return 0;

}
