#ifndef  PSPM_SOLVER_H_
#define  PSPM_SOLVER_H_

#include <vector>
#include "pspm_ode_solver2.h"

enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT};

template <class Model>
class Solver{
	private:
	Model * mod;				// Model should be a class with growth, mortality, birth, env, and IC functions
	
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
	
	public:
	double xb, xm;
	
	// The current state of the system, {t, S, dS/dt} 
	double current_time = 0;      // These are synced with ODE solver after each successful step
	std::vector <double> state;	  // +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 
		
	Solver(int _J, double _xb, double _xm, PSPM_SolverType _method);
	Solver(std::vector<double> xbreaks, PSPM_SolverType _method);

	void setModel(Model *M);
	void setInputNewbornDensity(double input_u0);
//	template<typename Func>
//	void initialize(Func calcIC);
	
	void initialize();

	const int size();
	const int xsize();
	const double* getX();
	vector<double> getx();

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

	double newborns_out();  // This is the actual system reproduction (fitness) hence biologically relevant
	double u0_out();        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value
	
	void print();
	
	template<typename wFunc>
	double integrate_x(wFunc w, double t, vector<double>&S, int power);
	
};

#include "../src/solver.tpp"

#endif



