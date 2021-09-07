#ifndef  PSPM_PSPM_SOLVER_H_
#define  PSPM_PSPM_SOLVER_H_

#include <vector>
#include <list>
#include <map>
#include <string>

#include "environment_base.h"
#include "species.h"
#include "pspm_ode_solver3.h"

enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT, SOLVER_IFMU};

class Solver{
	private:
	PSPM_SolverType method;
	RKCK45<vector<double>> odeStepper;

	int n_statevars_internal = 0;		// number of internal i-state variables (x and/or u)
	int n_statevars_system = 0;			// number of s-state variables 

	public:	
	EnvironmentBase * env;
	
	public:
	// The current state of the system, {t, S, dS/dt} 
	std::vector <double> state;		// +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 
	double current_time;			// these are synced with ODE solver only after a successful step

	public:
	std::vector<Species_Base*> species_vec;	

	struct{
		double ode_eps = 1e-6;	// FIXME: Need a function to set these params, so ODE solver can be reset
		double ode_initial_step_size = 1e-6;
		double convergence_eps = 1e-6;
		double cm_grad_dx = 1e-6;
		bool update_cohorts = true;
		int  max_cohorts = 500;
		double ebt_ucut = 1e-10;
		double ebt_grad_dx = 1e-6;
		std::string ode_method = "rk45ck";
		double ode_rk4_stepsize = 0.1;
		double ode_ifmu_stepsize = 0.1;
		bool ifmu_centered_grids = true;
	} control;
	
	bool use_log_densities = true;

	public:	
	Solver(PSPM_SolverType _method);

	void addSystemVariables(int _s);
	void addSpecies(int _J, double _xb, double _xm, bool log_breaks, Species_Base* _mod, int n_extra_vars, double input_birth_flux = -1);
	void addSpecies(std::vector<double> xbreaks, Species_Base* _mod, int n_extra_vars, double input_birth_flux = -1);

	//Species<Model>* get_species(int id);
	int n_species();

	void setEnvironment(EnvironmentBase * _env);

	//int setupLayout(Species<Model> &s);
	void resetState(); 	
	void resizeStateFromSpecies();

	void initialize();

	void copyStateToCohorts(std::vector<double>::iterator state_begin);		////const int size();
	void copyCohortsToState();
	
	////const int xsize();
	////const double* getX();
	////vector<double> getx();
	////double getMaxSize(vector<double>::iterator sbegin);
	double maxSize();
	//double get_u0(double t, int s);	


	////void calcRates_extra(double t, vector<double>&S, vector<double>& dSdt);
	
	void calcRates_FMU(double t, vector<double> &S, vector<double> &dSdt);
	
	void step_iFMU(double t, vector<double> &S, double dt);
	
	void calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt);
	void addCohort_EBT();
	void removeDeadCohorts_EBT();
	////vector<double> cohortsToDensity_EBT(vector <double> &breaks);

	void calcRates_CM(double t, vector<double>&S, vector<double> &dSdt);
	//double calc_u0_CM();
	void addCohort_CM();
	void removeCohort_CM();
	//void removeDenseCohorts_CM();
	
	
	void step_to(double tstop);

	////double stepToEquilibrium();

	void preComputeSpecies(int k, double t);
	double calcSpeciesBirthFlux(int k, double t);
	std::vector<double> newborns_out();  // This is the actual system reproduction (fitness) hence biologically relevant
	vector<double> u0_out();        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value
	////double get_u0_out();	// just returns from history without recomputing

	void print();
	
	// integrals over size and density
	template<typename wFunc>
	double integrate_x(wFunc w, double t, int species_id);

	template<typename wFunc>
	double integrate_wudx_above(wFunc w, double t, double xlow, int species_id);
	

	std::vector<double> getDensitySpecies_EBT(int k, int nbreaks);
};

//#include "../src/solver.tpp"
//#include "../src/mu.tpp"
//#include "../src/ebt.tpp"
//#include "../src/cm.tpp"
#include "../src/size_integrals.tpp"


#endif



