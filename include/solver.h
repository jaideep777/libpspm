#ifndef  PSPM_SOLVER_H_
#define  PSPM_SOLVER_H_

#include <vector>
#include <list>
#include <map>

#include "pspm_ode_solver2.h"
#include "iterator_set.h"

enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT};

template <class Model>
class Solver{
	private:
	Model * mod;				// Model should be a class with growth, mortality, birth, env, and IC functions
	
	std::vector<string> varnames;				// state has internal variables (x, u) and possibly extra variables 
	std::vector<string> varnames_extra;   //   +-- which will be created in the state in this order (for each species)
	std::vector<int> strides;				// defines how the variables are packed into the state vector
	std::vector<int> offsets;				// 


	PSPM_SolverType method;

	int J;
	int state_size;
	//int nx;
	
	// These are only initial values, and repeatedly used only by FMU. 
	// Others derive them from the state. 
	// Kept private so users dont accidently access them for other solvers
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <double> h;
	
	RKCK45<vector<double>> odeStepper;

	double u0_in = -1;
	
	std::list<double> u0_out_history;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	std::vector<double> schedule;

	public:
	double xb, xm;
	
	// The current state of the system, {t, S, dS/dt} 
	double current_time;			// These are synced with ODE solver after each successful step
	std::vector <double> state;		// +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 

	public:	
	Solver(int _J, double _xb, double _xm, PSPM_SolverType _method);
	Solver(std::vector<double> xbreaks, PSPM_SolverType _method);
	void resetState(const std::vector<double>& xbreaks); 	

	void setModel(Model *M);
	void setInputNewbornDensity(double input_u0);
	double createSizeStructuredVariables(std::vector<std::string> names);
	
	void initialize();

	const int size();
	const int xsize();
	const double* getX();
	vector<double> getx();
	IteratorSet<vector<double>::iterator> getIterators_state();
	IteratorSet<vector<double>::iterator> getIterators_rates();

	void calcRates_extra(double t, vector<double>&S, vector<double>& dSdt);
	
	void calcRates_FMU(double t, vector<double> &U, vector<double> &dUdt);
	//void calcRates_FMU(double t);	
	
	void calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt);
	void addCohort_EBT();
	void removeDeadCohorts_EBT();
	vector<double> cohortsToDensity_EBT(vector <double> &breaks);

	void calcRates_CM(double t, vector<double>&S, vector<double> &dSdt);
	double calc_u0_CM();
	void addCohort_CM();
	void removeCohort_CM();
	
	
	void step_to(double tstop);

	double stepToEquilibrium();

	double newborns_out();  // This is the actual system reproduction (fitness) hence biologically relevant
	double u0_out();        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value
	double get_u0_out();	// just returns from history without recomputing

	void print();
	
	template<typename wFunc>
	double integrate_x(wFunc w, double t, vector<double>&S, int power);
	
};

#include "../src/solver.tpp"

#endif



