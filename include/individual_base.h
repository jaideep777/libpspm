#ifndef  PSPM_PSPM_INDIVIDUAL_BASE_H_
#define  PSPM_PSPM_INDIVIDUAL_BASE_H_

#include <iostream>
#include <vector>
#include <string>

class IndividualBase{
	public:
	std::vector<std::string> varnames;
	
	// essential functions which must be defined by user
	virtual double init_density(double x, void * _env, double bf) = 0;
	virtual double growthRate(double x, double t, void * _env) = 0;
	virtual double mortalityRate(double x, double t, void * _env) = 0;
	virtual double birthRate(double x, double t, void * _env) = 0;

	virtual ~IndividualBase();
	
	virtual void   set_size(double _x);
	virtual void   preCompute(double x, double t, void * _env);
	virtual double establishmentProbability(double t, void  * _env);

	// Functions related to additional state variables
	virtual void init_state(double t, void * _env);
	virtual std::vector<double>::iterator set_state(std::vector<double>::iterator &it);
	virtual std::vector<double>::iterator get_state(std::vector<double>::iterator &it);
	virtual std::vector<double>::iterator get_rates(std::vector<double>::iterator &it); 

	// printing
	virtual void print(std::ostream &out = std::cout);
};


#endif


