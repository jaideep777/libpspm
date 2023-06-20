#ifndef  PSPM_PSPM_SPECIES_H_
#define  PSPM_PSPM_SPECIES_H_

#include <vector>
#include <list>
#include <string>

#include "tensor.h"
#include "tensor_cohort.h"

//#include "iterator_set.h"
#include "cohort.h"

// forward declaration of Solver so Species can befriend it
class Solver;


class Species_Base{
	// Solver should be able to access Species' privates
	friend class Solver;

	protected: // private members
	int J;	
	int n_extra_statevars = 0;

	std::list<double> birth_flux_out_history;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	
	// These are only used by FMU solver. 
	// Others derive them from the state. 
	// Kept private so users dont accidently access them for other solvers
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <std::vector<double>> Xn; // multiple state model
	std::vector <std::vector<double>> xn; // multiple state model
	std::vector <double> h;
	// std::vector <double> schedule; // used only by CM/EBT

	double noff_abm = 0; // used by ABM solver to insert offspring

	public:
	double birth_flux_in;
	
	//debug only
	bool bfin_is_u0in = false;
	
	double xb; 
	std::vector <double> xnb;
	
	public: // public functions
	// virtual Species_Base * create() = 0;
	virtual ~Species_Base() = 0;
	
	int xsize(int k);
	int cohortsize();
	int statesize();

	int xsize();
	int size();

	void set_inputBirthFlux(double b);
	void set_bfin_is_u0in(bool flag);

	public:
	virtual void resize(int _J) = 0;
	virtual double get_maxSize() = 0;
	virtual void print() = 0;
	virtual void print_extra(); // not pure virtual, by defualt, there is nothing extra to print.

	virtual void set_xb(double _xb) = 0;
	virtual void set_xnb(std::vector<double> _xnb) = 0;
	virtual void set_ub(double _ub) = 0;
	virtual void set_birthTime(int i, double t0) = 0;
	virtual void setX(int i, double _x) = 0;
	virtual void setXn(int i, int k, double _x) = 0;
	virtual void setXn(int i, std::vector<double> _x) = 0;
	virtual void setU(int i, double _u) = 0;

	virtual double getX(int i) = 0;
	virtual std::vector<double> getXn(int i) = 0;
	virtual double getU(int i) = 0;

	virtual std::vector <double> getStateAt(int i) = 0;
	virtual double dXn (int i) = 0;

	virtual double init_density(int i, double x, void * env) = 0;
	virtual double init_density(int i, std::vector<double> x, void * env) = 0;
	virtual void initExtraState(double t, void * env) = 0;
	virtual void initAndCopyExtraState(double t, void * env, std::vector<double>::iterator &it) = 0;
	virtual void initBoundaryCohort(double t, void * env) = 0;

	virtual void copyExtraStateToCohorts(std::vector<double>::iterator &it) = 0;
	virtual void copyCohortsExtraToState(std::vector<double>::iterator &it) = 0;
	
	virtual double establishmentProbability(double t, void * env) = 0;
	virtual double calc_boundary_u(double gb, double pe) = 0;
	virtual double get_boundary_u() = 0;

	virtual void triggerPreCompute() = 0;

	// TODO: argument x can probably be removed from these functions
	virtual double growthRate(int i, double x, double t, void * env) = 0;
	virtual double growthRateOffset(int i, double x, double t, void * env) = 0;
	virtual std::vector<double> growthRateGradient(int i, double x, double t, void * env, double grad_dx) = 0;
	virtual std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env) = 0;
	virtual double mortalityRate(int i, double x, double t, void * env) = 0;
	virtual std::vector<double> mortalityRateGradient(int i, double x, double t, void * env, double grad_dx) = 0;
	virtual double birthRate(int i, double x, double t, void * env) = 0;
	virtual void getExtraRates(std::vector<double>::iterator &it) = 0;

	virtual void addCohort(int n = 1) = 0;
	template<class T> void addCohort(T bc);

	virtual void markCohortForRemoval(int i) = 0;
	virtual void removeMarkedCohorts() = 0;
	virtual void removeDensestCohort() = 0;
	virtual void removeDenseCohorts(double dxcut) = 0;
	virtual void removeDeadCohorts(double ucut) = 0;
	virtual void mergeCohortsAddU(double dxcut) = 0;

	virtual void sortCohortsDescending(int skip=0) = 0;
	
	virtual void save(std::ofstream &fout) = 0;
	virtual void restore(std::ifstream &fin) = 0;

//	virtual void backupCohort(int j) = 0;
//	virtual void restoreCohort(int j) = 0;
//	virtual void copyBoundaryCohortTo(int j) = 0;
};



template <class Model>
class Species : public Species_Base{
	protected:
	TensorCohort<Cohort<Model>> cohorts;
	Cohort<Model> boundaryCohort;
	
	//Cohort<Model> savedCohort; // a cohort to save a backup of any other cohort

	public:
	// TODO: make these virtual? - not needed. They are virtual by default.
	Species(std::vector<double> breaks = std::vector<double>());
	Species(Model M);
	// Species<Model>* create(); // virtual constructor needed for deserialization

	
	void resize(int _J);
	void resize(std::vector<int> xsize);
	double get_maxSize();
	
	void print();
	
	void set_xb(double _xb);
	void set_xnb(std::vector <double> _xnb);
	void set_ub(double _ub);
	void set_birthTime(int i, double t0);
	void setX(int i, double _x);
	void setXn(int i, int k, double _x);
	void setXn(int i, std::vector<double> _x);
	void setU(int i, double _u);
	
	double getX(int i);
	std::vector<double> getXn(int i);
	double getU(int i);

	std::vector <double> getStateAt(int i);
	double dXn (int i);
	
	double init_density(int i, double x, void * env);
	void initExtraState(double t, void * env);
	void initAndCopyExtraState(double t, void * env, std::vector<double>::iterator &it);
	void initBoundaryCohort(double t, void * env);

	void copyExtraStateToCohorts(std::vector<double>::iterator &it);
	void copyCohortsExtraToState(std::vector<double>::iterator &it);

	double establishmentProbability(double t, void * env);
	double calc_boundary_u(double gb, double pe);
	double get_boundary_u();
	
	void triggerPreCompute();
	double growthRate(int i, double x, double t, void * env);
	double growthRateOffset(int i, double x, double t, void * env);
	std::vector<double> growthRateGradient(int i, double x, double t, void * env, double grad_dx);
	std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env);
	double mortalityRate(int i, double x, double t, void * env);
	std::vector<double> mortalityRateGradient(int i, double x, double t, void * env, double grad_dx);
	double birthRate(int i, double x, double t, void * env);
	void getExtraRates(std::vector<double>::iterator &it);

	void addCohort(int n = 1);
	void addCohort(Cohort<Model> bc);

	void markCohortForRemoval(int i);
	void removeMarkedCohorts();
	void removeDensestCohort();
	void removeDenseCohorts(double dxcut);
	void removeDeadCohorts(double ucut);
	void mergeCohortsAddU(double dxcut);

	void sortCohortsDescending(int skip=0);
	
	void save(std::ofstream &fout);
	void restore(std::ifstream &fin);

//	void backupCohort(int j);
//	void restoreCohort(int j);
//	void copyBoundaryCohortTo(int j);

	public:
	Cohort<Model>& getCohort(int i);
};


#include "../src/species.tpp"

#endif
