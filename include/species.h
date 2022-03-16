#ifndef  PSPM_PSPM_SPECIES_H_
#define  PSPM_PSPM_SPECIES_H_

#include <vector>
#include <list>
#include <string>

//#include "iterator_set.h"
#include "cohort.h"

// forward declaration of Solver so Species can befriend it
class Solver;


class Species_Base{
	// All kinds of Solvers should be friends of Species
	friend class Solver;

	protected: // private members
	//int start_index;
	int J;	
	int n_extra_statevars = 0;
	//std::vector<std::string> varnames;			// state has internal variables (x, u) and possibly extra variables 
	//std::vector<std::string> varnames_extra;		//   +-- which will be created in the state in this order (for each species)
	//std::vector<int> strides;				// defines how the variables are packed into the state vector
	//std::vector<int> offsets;				// 
	

	std::list<double> birth_flux_out_history;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	
	// These are only used by FMU solver. 
	// Others derive them from the state. 
	// Kept private so users dont accidently access them for other solvers
	std::vector <double> X;	
	std::vector <double> x;
	std::vector <double> h;
	std::vector <double> schedule; // used only by CM/EBT

	//public:
	//double u0_save;

	public:
	double birth_flux_in;
	
	//debug only
	bool bfin_is_u0in = false;
	
	//private: // private functions
	//int addVar(std::string name, int stride, int offset);
	//void clearVars();

	public: // public members
	double xb; //, xm;  // FIXME. Dangerous because xm is never updated
	bool is_resident;
	
	//Model * mod = nullptr;
	public:

	public: // public functions

	virtual ~Species_Base() = 0;
	
	int xsize();
	int size();

	void set_inputBirthFlux(double b);
	void set_bfin_is_u0in(bool flag);

	public:
	virtual void resize(int _J) = 0;
	virtual double get_maxSize() = 0;
	virtual void print() = 0;

	virtual void set_xb(double _xb) = 0;
	virtual void set_birthTime(int i, double t0) = 0;
	virtual void setX(int i, double _x) = 0;
	virtual void setU(int i, double _u) = 0;

	virtual double getX(int i) = 0;
	virtual double getU(int i) = 0;

	virtual double init_density(int i, double x, void * env) = 0;
	virtual void initAndCopyExtraState(double t, void * env, std::vector<double>::iterator &it) = 0;
	virtual void initBoundaryCohort(double t, void * env) = 0;

	virtual void copyExtraStateToCohorts(std::vector<double>::iterator &it) = 0;
	virtual void copyCohortsExtraToState(std::vector<double>::iterator &it) = 0;
	
	virtual double establishmentProbability(double t, void * env) = 0;
	virtual double get_u0(double t, void * env) = 0;
	virtual double get_boundary_u() = 0;

	// TODO: argument x can probably be removed from these functions
	virtual void triggerPreCompute() = 0;
	virtual void preComputeAllCohorts(double t, void * env) = 0;
	virtual void preCompute(int i, double t, void * env) = 0;
	virtual double growthRate(int i, double x, double t, void * env) = 0;
	virtual double growthRateOffset(int i, double x, double t, void * env) = 0;
	virtual std::vector<double> growthRateGradient(int i, double x, double t, void * env, double grad_dx) = 0;
	virtual std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env) = 0;
	virtual double mortalityRate(int i, double x, double t, void * env) = 0;
	virtual std::vector<double> mortalityRateGradient(int i, double x, double t, void * env, double grad_dx) = 0;
	virtual double birthRate(int i, double x, double t, void * env) = 0;
	virtual void getExtraRates(std::vector<double>::iterator &it) = 0;

	virtual void addCohort() = 0;
	template<class T> void addCohort(T bc);

	virtual void removeDensestCohort() = 0;
	virtual void removeDenseCohorts(double dxcut) = 0;
	virtual void removeDeadCohorts(double ucut) = 0;

	virtual void backupCohort(int j) = 0;
	virtual void restoreCohort(int j) = 0;
	virtual void copyBoundaryCohortTo(int j) = 0;
};



template <class Model>
class Species : public Species_Base{
	protected:
	std::vector<Cohort<Model>> cohorts;
	Cohort<Model> boundaryCohort;
	
	Cohort<Model> savedCohort; // a cohort to save a backup of any other cohort

	public:
	// TODO: make these virtual?
	Species(std::vector<double> breaks = std::vector<double>());
	Species(Model M);
	
	void resize(int _J);
	double get_maxSize();
	
	void print();
	
	void set_xb(double _xb);
	void set_birthTime(int i, double t0);
	void setX(int i, double _x);
	void setU(int i, double _u);
	
	double getX(int i);
	double getU(int i);
	
	double init_density(int i, double x, void * env);
	void initAndCopyExtraState(double t, void * env, std::vector<double>::iterator &it);
	void initBoundaryCohort(double t, void * env);

	void copyExtraStateToCohorts(std::vector<double>::iterator &it);
	void copyCohortsExtraToState(std::vector<double>::iterator &it);

	double establishmentProbability(double t, void * env);
	double get_u0(double t, void * env);
	double get_boundary_u();
	
	void triggerPreCompute();
	void preComputeAllCohorts(double t, void * env);
	void preCompute(int i, double t, void * env);
	double growthRate(int i, double x, double t, void * env);
	double growthRateOffset(int i, double x, double t, void * env);
	std::vector<double> growthRateGradient(int i, double x, double t, void * env, double grad_dx);
	std::vector<double> growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env);
	double mortalityRate(int i, double x, double t, void * env);
	std::vector<double> mortalityRateGradient(int i, double x, double t, void * env, double grad_dx);
	double birthRate(int i, double x, double t, void * env);
	void getExtraRates(std::vector<double>::iterator &it);

	void addCohort();
	void addCohort(Cohort<Model> bc);

	void removeDensestCohort();
	void removeDenseCohorts(double dxcut);
	void removeDeadCohorts(double ucut);
	
	void backupCohort(int j);
	void restoreCohort(int j);
	void copyBoundaryCohortTo(int j);

	public:
	Cohort<Model>& getCohort(int i);
};


#include "../src/species.tpp"

#endif
