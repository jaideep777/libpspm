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
	
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <double> h;
	
	double newborns;
		
	RKCK45<vector<double>> odeStepper;

	public:
	double xb, xm;
	
	double current_time = 0;
	std::vector <double> state;
	std::vector <double> rates; 
		
	Solver(int _J, double _xb, double _xm, PSPM_SolverType _method);
	Solver(std::vector<double> xbreaks, PSPM_SolverType _method);

	void setModel(Model *M);
//	template<typename Func>
//	void initialize(Func calcIC);
	
	void initialize();

	const int size();
	const int xsize();
	const double* getX();
	vector<double> getx();

	void calcRates_FMU(double t, vector<double> &U, vector<double> &dUdt);
	//void calcRates_FMU(double t);	
	
	void addCohort_EBT();
	void removeDeadCohorts_EBT();
	void calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt);
	vector<double> cohortsToDensity_EBT(vector <double> &breaks);

	void calcRates_CM(double t, vector<double>&S, vector<double> &dSdt);
	double calcBirthFlux_CM(double _u0);
	void addCohort_CM(double u0 = -1);
	void removeCohort_CM();
	
	
	void step_to(double tstop);

	void print();
	
	template<typename wFunc>
	double integrate_x(wFunc w, int power);
	
};

#include "../src/solver.tpp"

#endif



