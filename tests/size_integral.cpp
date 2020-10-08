#include <solver.h>
#include "test_model.h"

int main(){
	Solver<TestModel> S({0,1,2,3}, SOLVER_CM);
	//S.print();

	S.state = {0,1,2,3, log(1),log(2),log(4),log(8)};
	//S.print();

	auto w = [](double x, double t){return 1;};

	cout << S.integrate_x(w, 0, S.state, 1) << " " << S.integrate_wudx_above(w, 0, 0, S.state) << endl;
	
	if(abs(S.integrate_x(w, 0, S.state, 1) - 10.5) > 1e-5) return 1;
	if(abs(S.integrate_x(w, 0, S.state, 1) - S.integrate_wudx_above(w, 0, 0, S.state)) > 1e-5) return 1;
	//if(abs(S.integrate_wudx_above(w, 0, 1.5, S.state) - 7.75) < 2e-5) return 1;

	return 0;

}
