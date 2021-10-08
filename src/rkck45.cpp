#include "rkck45.h"
using namespace std;


RKCK45::RKCK45(double t_start_, double accuracy, double h1) :
		ht(h1), eps_rel(accuracy), eps_abs(accuracy), xt(t_start_){
		
}

void RKCK45::resize(int new_size){
	sys_size = new_size;

	yscal.resize(new_size);
	dydx.resize(new_size);
	k1.resize(new_size);
	k2.resize(new_size);
	k3.resize(new_size);
	k4.resize(new_size);
	k5.resize(new_size);
	yt.resize(new_size);
}

