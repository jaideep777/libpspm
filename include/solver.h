#ifndef  PSPM_SOLVER_H_
#define  PSPM_SOLVER_H_

#include <vector>


enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT};

template <class Model>
class Solver{
	private:
	Model * mod;				// Model should be a class with growth, mortality, birth, env, and IC functions
	
	PSPM_SolverType method;

	int J;
	int state_size;
	int nx;
	
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <double> h;

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

	void calcRates_FMU(double t);	
	void calcRates_EBT();	
	
	void print();
	
	template<typename wFunc>
	double integrate_x(wFunc w, int power);
	
	void step_to(double tf);
};

#include "../src/solver.tpp"

#endif



