#include <solver.h>
#include "test_model.h"

int main(){
	Solver<TestModel> S({0,1,2,3}, SOLVER_CM);
	//S.print();

	S.state = {0,1,2,3, 1,2,4,8};
	//S.print();

	auto w = [](double x, double t){return 1;};

	assert(S.integrate_x(w, 0, S.state, 1) == 10.5);
	assert(S.integrate_x(w, 0, S.state, 1) == S.integrate_wudx_above(w, 0, 0, S.state));
	assert(S.integrate_wudx_above(w, 0, 1.5, S.state) == 7.75);

	return 0;

}
