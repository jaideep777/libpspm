#include "rkck45.h"
#include "io_utils.h"
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


void RKCK45::save(std::ofstream &fout){
	fout << "OdeSolver::v1\n";
	fout << std::make_tuple(
			ht
			, eps_rel
			, eps_abs 
			, xt
			, t_stop
			, nok
			, nbad
			, nfe
			, sys_size);
	fout << '\n';
	fout << yscal;
}

void RKCK45::restore(std::ifstream &fin){
	string s; fin >> s; // version number (discard)
	fin >> ht
	    >> eps_rel
	    >> eps_abs 
	    >> xt
	    >> t_stop
	    >> nok
	    >> nbad
	    >> nfe
	    >> sys_size;
	resize(sys_size);
	fin >> yscal;
}

