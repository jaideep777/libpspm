#ifndef  PSPM_PSPM_INDIVIDUAL_BASE_H_
#define  PSPM_PSPM_INDIVIDUAL_BASE_H_

#include <iostream>
#include <vector>

class Individual_Base{
	public:
	virtual ~Individual_Base() = 0;
	virtual double init_density(double x, void * _env, double bf) = 0;
	virtual void preCompute(double x, double t, void * _env) = 0;
	virtual double growthRate(double x, double t, void * _env) = 0;
	virtual double mortalityRate(double x, double t, void * _env) = 0;
	virtual double birthRate(double x, double t, void * _env) = 0;
	virtual double establishmentProbability(double t, void  * _env) = 0;
	virtual void set_size(double _x) = 0;
	virtual void init_state(double t, void * env) = 0;
	virtual std::vector<double>::iterator set_state(std::vector<double>::iterator &it) = 0;
	virtual std::vector<double>::iterator get_state(std::vector<double>::iterator &it) = 0;
	virtual std::vector<double>::iterator get_rates(std::vector<double>::iterator &it) = 0; 
	virtual void print(std::ostream &out = std::cout) = 0;
};

inline Individual_Base::~Individual_Base(){}

#endif


