#ifndef  PSPM_SOLVER_H_
#define  PSPM_SOLVER_H_

#include <vector>


enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT};

class Solver{
	private:
	PSPM_SolverType method;
	int J;
	int state_size;
	int nx;
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <double> h;

	public:
	double xb, xm;

	std::vector <double> state;
	std::vector <double> rates; 
		
	Solver(int _J, double _xb, double _xm, PSPM_SolverType _method);
	Solver(std::vector<double> xbreaks, PSPM_SolverType _method);

//	template<typename Func>
//	void initialize(Func calcIC);
	
	template<class _Model>
	void initialize(_Model mod);

	const int size();
	const int xsize();
	const double* getX();

//	template<typename growthFunc, typename mortFunc, typename fecFunc, typename envFunc>
//	void calcRates_FMU(); // (Func w)	
	
	
	void print();
};

#include "../src/solver.tpp"

#endif



