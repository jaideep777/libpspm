#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "gsl_ode_solver.h"
#include "solver.h"

using namespace std;

int main(){
	
	Solver S(10, 0,1, SOLVER_FMU);
	S.print();

	return 0;
}


int main(){

	Environment env;
	
	Solver S(25, "FMU");

	int rate_eval_count = 0;
	auto func_lambda = [&S, &env, &rate_eval_count](double t, const double y[], double f[]) -> int {
		env.time = t;
		S.setState(y);
		env.value = computeEnvironment(S.X, S.U, S.h);
//		cout << env.value << endl;
		S.getRates(env, f);

		++rate_eval_count;
		return GSL_SUCCESS;
	};


	ODE_Solver<decltype(func_lambda)> sys(func_lambda, S.J);

	vector <double> state0 = initialize(S.X);
	cout << "initial state: "; 
	for (auto x : state0) cout << x << " ";
	cout << endl;
	
	sys.initialize(state0, 0);
	S.U = state0;

	double t0 = 0, tf = 5, dt = .1;
	size_t nsteps = (tf-t0)/dt;


	ofstream fout("patch_full_hts.txt");

	cout << "Step " << 0 <<  endl;
	fout << env.time << "\t";
	for (int i=0; i<S.J; ++i){		
		fout << S.U[i] << "\t";// <<
	}						
	fout << endl;

	for (size_t i=1; i < nsteps; ++i){
		
		cout << "Step " << i << endl;
		sys.step_to(t0+i*dt);

		fout << t0+i*dt << "\t";
		for (int j=0; j<S.J; ++j){		
			fout << S.U[j] << "\t";// <<
		}						
		fout << endl;

	}

	fout.close();

	return 0;
}


