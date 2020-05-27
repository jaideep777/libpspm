#include <vector>
#include <iostream>
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
	int m;
	cout << "Enter method: ";
	cin >> m;

	Solver S(10, 0,1, PSPM_SolverType(m));
	
	double sc = 2;
	auto f = [sc](double x){
		return 10*x/sc;
	};

	Model M;

	S.initialize(M);
	S.print();	

	
	return 0;
}


