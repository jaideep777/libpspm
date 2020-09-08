#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

class Model{
	public:
	double sc = 10;
	
	double calcIC(double x){
		return sc*x;
	}
};


int main(){

	Solver<Model> S(10, 0,1, SOLVER_FMU);
	
	Model M;
	
	S.setModel(&M);
	S.initialize();
	//S.print();	
	
	for (int i=0; i<10; ++i){
		float X = .1*i+0.05;
		if (fabs(S.state[i] - M.calcIC(X)) > 1e-6) return 1;
	}
	
	return 0;
}


