#ifndef  PSPM_PSPM_SOLVER_H_
#define  PSPM_PSPM_SOLVER_H_

#include <vector>
#include <list>
#include <map>
#include <string>
#include <random>

#include "environment_base.h"
#include "species.h"
#include "ode_solver.h"
#include "cubic_spline.h"

enum PSPM_SolverType {SOLVER_FMU, SOLVER_MMU, SOLVER_CM, SOLVER_EBT, SOLVER_IFMU, SOLVER_ABM, SOLVER_IEBT};

class Solver{
	private:
	static std::map<std::string, PSPM_SolverType> methods_map;

	static constexpr bool debug = false;
	PSPM_SolverType method;

	int n_statevars_internal = 0;		// number of internal i-state variables (x and/or u)
	int n_statevars_system = 0;			// number of s-state variables 

	std::default_random_engine generator; // random number generator

	public:	
	OdeSolver odeStepper;
	EnvironmentBase * env = nullptr;
	
	// The current state of the system, {t, S, dS/dt} 
	double current_time = -999;			// these are synced with ODE solver only after a successful step
	std::vector <double> state;		// +-- They are NOT synced during the ODE solver's internal steps
	std::vector <double> rates; 

	std::vector<Species_Base*> species_vec;	

	struct{
		double ode_eps = 1e-6;	// FIXME: Need a function to set these params, so ODE solver can be reset
		double ode_initial_step_size = 1e-6;
		double convergence_eps = 1e-6;
		double cm_grad_dx = 1e-6;
		bool update_cohorts = true;
		bool cm_remove_cohorts = true;
		int  max_cohorts = 500;
		double cm_dxcut = 1e-10;
		double ebt_ucut = 1e-10;
		double ebt_grad_dx = 1e-6;
		double ode_rk4_stepsize = 0.1;
		double ode_ifmu_stepsize = 0.1;
		bool ifmu_centered_grids = true;
		bool integral_interpolate = true;
		double ifmu_order = 1;
		double abm_n0 = 100;
		double abm_stepsize = 0.02;
	} control;
	
	bool use_log_densities = true;

	private:
	double pi0, N0;
	void realizeEbtBoundaryCohort(Species_Base * spp);
	void restoreEbtBoundaryCohort(Species_Base * spp);

	public:	
	Solver(PSPM_SolverType _method, std::string ode_method = "rk45ck");
	Solver(std::string _method, std::string ode_method = "rk45ck");

	void addSystemVariables(int _s);
	void addSpecies(int _J, double _xb, double _xm, bool log_breaks, Species_Base* _mod, int n_extra_vars, double input_birth_flux = -1);
	void addSpecies(std::vector<double> xbreaks, Species_Base* _mod, int n_extra_vars, double input_birth_flux = -1);
	void removeSpecies(Species_Base* spp);

	//Species<Model>* get_species(int id);
	int n_species();

	void setEnvironment(EnvironmentBase * _env);

	void resetState(double t0 = 0);
	void resizeStateFromSpecies();

	void initializeSpecies(Species_Base * s);
	void initialize();

	void copyStateToCohorts(std::vector<double>::iterator state_begin);		////const int size();
	void copyCohortsToState();
	
	double maxSize();

	void updateEnv(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);

	/// @brief calculate \f$du/dt\f$ using the FMU solver
	void calcRates_FMU(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	
	void calcOdeRatesImplicit(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);

	void stepU_iFMU(double t, std::vector<double> &S, std::vector<double> &dSdt, double dt);
	void stepU_iEBT(double t, std::vector<double> &S, std::vector<double> &dSdt, double dt);

	/// @brief calculate \f$du/dt\f$ using the EBT solver
	void calcRates_EBT(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	void addCohort_EBT();
	void removeDeadCohorts_EBT();

	/// @brief calculate \f$du/dt\f$ using the CM solver
	void calcRates_CM(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt);
	void addCohort_CM();
	void removeCohort_CM();

	void stepABM(double t, double dt);
	
	template<typename AfterStepFunc>
	void step_to(double tstop, AfterStepFunc &afterStep_user);

	void step_to(double tstop);

	double calcSpeciesBirthFlux(int k, double t); // TODO: Make this private. Because this does not update Env, so unsafe for users to call. Users should call newborns_out() instead because it does updateEnv + calcSpeciesBirthFlux()
	std::vector<double> newborns_out(double t);  // This is the actual system reproduction (fitness) hence biologically relevant
	std::vector<double> u0_out(double t);        // This is used for equilibrium solving, because in general, u0 rather than birthFlux, will approach a constant value

	void print();
	
	// integrals over size and density
	/// @brief Computes the integral \f[I = \int_{x_b}^{x_m} w(z,t)u(z)dz\f] for the specified species. For details, see @ref integrate_wudx_above
	template<typename wFunc>
	double integrate_x(wFunc w, double t, int species_id);

	/// @brief Computes the partial integral \f[I = \int_{x_{low}}^{x_m} w(z,t)u(z)dz\f] for the specified species. 
	template<typename wFunc>
	double integrate_wudx_above(wFunc w, double t, double xlow, int species_id);

	std::vector<double> getDensitySpecies(int k, std::vector<double> breaks, Spline::Extr extrapolation_method = Spline::ZERO);

	void save(std::ofstream &fout);
	void restore(std::ifstream &fin, std::vector<Species_Base*> spp_proto);

};

#include "../src/solver.tpp"
#include "../src/size_integrals.tpp"


#endif



