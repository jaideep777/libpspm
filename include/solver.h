#ifndef  PSPM_PSPM_SOLVER_H_
#define  PSPM_PSPM_SOLVER_H_

#include <vector>
#include <list>
#include <map>

#include "species.h"
#include "pspm_ode_solver3.h"
#include "iterator_set.h"

enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT};

template <class Model, class Environment>
class Solver{
	private:
	PSPM_SolverType method;
	RKCK45<vector<double>> odeStepper;
	
	std::vector<Species<Model>> species_vec;	
	Environment * env;

	public:
	// The current state of the system, {t, S, dS/dt} 
	double current_time;			// These are synced with ODE solver after each successful step
	std::vector <double> state;		// +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 

	struct{
		double ode_eps = 1e-6;	// FIXME: Need a function to set these params, so ODE solver can be reset
		double ode_initial_step_size = 1e-6;
		double convergence_eps = 1e-6;
		double cm_grad_dx = 1e-6;
	} control;
	
	bool use_log_densities = true;

	public:	
	Solver(PSPM_SolverType _method);

	void addSpecies(int _J, double _xb, double _xm, bool log_breaks, Model* _mod, std::vector<std::string> extra_vars = std::vector<std::string>(), double input_birth_flux = -1);
	void addSpecies(std::vector<double> xbreaks, Model* _mod, std::vector<std::string> extra_vars = std::vector<std::string>(), double input_birth_flux = -1);

	Species<Model>* get_species(int id);
	int n_species();

	void setEnvironment(Environment * _env);

	int setupLayout(Species<Model> &s);
	void resetState(); 	

	void initialize();

	//const int size();
	//const int xsize();
	//const double* getX();
	//vector<double> getx();
	//double getMaxSize(vector<double>::iterator sbegin);
	double maxSize(std::vector<double>::iterator state_begin);
	double get_u0(double t, int s);	


	//void calcRates_extra(double t, vector<double>&S, vector<double>& dSdt);
	
	void calcRates_FMU(double t, vector<double> &S, vector<double> &dSdt);
	
	void calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt);
	void addCohort_EBT();
	void removeDeadCohorts_EBT();
	//vector<double> cohortsToDensity_EBT(vector <double> &breaks);

	void calcRates_CM(double t, vector<double>&S, vector<double> &dSdt);
	double calc_u0_CM();
	void addCohort_CM();
	void removeCohort_CM();
	
	
	void step_to(double tstop);

	//double stepToEquilibrium();

	std::vector<double> newborns_out();  // This is the actual system reproduction (fitness) hence biologically relevant
	vector<double> u0_out();        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value
	//double get_u0_out();	// just returns from history without recomputing

	void print();
	
	// integrals over size and density
	template<typename wFunc>
	double integrate_x(wFunc w, double t, vector<double>&S, int species_id);

	template<typename wFunc>
	double integrate_wudx_above(wFunc w, double t, double xlow, vector<double>&S, int species_id);
	
};

#include "../src/solver.tpp"
#include "../src/mu.tpp"
#include "../src/ebt.tpp"
#include "../src/cm.tpp"
#include "../src/size_integrals.tpp"


#endif



