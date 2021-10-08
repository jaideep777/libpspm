#ifndef PSPM_ODE_SOLVER_H_
#define PSPM_ODE_SOLVER_H_

#include "rkck45.h"
#include "lsoda.h"

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

enum SolverType {ODE_RKCK45, ODE_LSODA};

class OdeSolver{
	private:
	SolverType type;
	void * solver = nullptr;
	int nfe_cumm = 0;

	public:
	struct{
		double abs_tol = 1e-8;
		double rel_tol = 1e-8;
	} control;

	public:
	OdeSolver(std::string method, double t_start, double rtol, double atol){
		if (method == "rk45ck") type = ODE_RKCK45;
		else if (method == "lsoda") type = ODE_LSODA;
		else {
			std::cout << "Fatal: Unknown ODE method " << method << "\n";
			exit(1);
		}
		reset(t_start, rtol, atol);
	}

	~OdeSolver(){
		if      (type == ODE_RKCK45) delete static_cast<RKCK45*>(solver);	
		else if (type == ODE_LSODA)  delete static_cast<LSODA*>(solver);	
	}

	
	void reset(double t_start, double rtol, double atol){
		nfe_cumm = 0;
		control.abs_tol = atol;
		control.rel_tol = rtol;
		if (type == ODE_RKCK45){
			RKCK45* sol = static_cast<RKCK45*>(solver);
			delete sol; sol = nullptr;
			solver = new RKCK45(t_start, rtol, 1e-8);
		}
		else if (type == ODE_LSODA){
			LSODA* sol = static_cast<LSODA*>(solver);
			delete sol; sol = nullptr;
			solver = new LSODA();
		}
	}


	template <class Functor>
	void step_to(double t_stop, double &t, std::vector<double>&y, Functor &derivs){
		if (t_stop == t || y.size() == 0) return;
		
		if (type == ODE_RKCK45){
			(static_cast<RKCK45*>(solver))->Step_to(t_stop, t, y, derivs);
		}
		else if (type == ODE_LSODA ){
			LSODA* sol = static_cast<LSODA*>(solver);
			sol->set_istate(1); // forces re-initialization
			sol->lsoda_update(derivs, y.size(), y, &t, t_stop, nullptr, control.rel_tol, control.abs_tol);
			if(sol->get_istate() <= 0) {
				std::cerr << "LSODA Error: istate = " << sol->get_istate() << std::endl;
				return;
			}
			nfe_cumm += sol->get_fncalls();
		}
		
		//else if (type == ODE_LSODA ){
			////LSODA* sol = (LSODA*)solver;
			////delete sol;
			//LSODA * sol = new LSODA();
			//sol->lsoda_update(derivs, y.size(), y, &t, t_stop, nullptr, control.rel_tol, control.abs_tol);
			//if(sol->get_istate() <= 0) {
				//std::cerr << "LSODA Error: istate = " << sol->get_istate() << std::endl;
				//return;
			//}
			//nfe_cumm += sol->get_fncalls();
			//delete sol;
		//}
	}

	int get_fn_evals(){
		if      (type == ODE_RKCK45) return static_cast<RKCK45*>(solver)->get_fn_evals();	
		else if (type == ODE_LSODA)  return nfe_cumm;	
	}

};



#endif

