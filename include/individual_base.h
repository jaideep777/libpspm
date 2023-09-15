#ifndef  PSPM_PSPM_INDIVIDUAL_BASE_H_
#define  PSPM_PSPM_INDIVIDUAL_BASE_H_

#include <iostream>
#include <vector>
#include <array>
#include <string>

#include "io_utils.h"

template<size_t dim = 1>
class IndividualBase{
	public:

	// Since we cant force users to manage x, these variables are managed by Cohort.
	// But it is defined here becuase it is templated and Cohort doesnt have access to `dim`
	size_t state_size = dim;
	std::array<double, dim> x; 
	// ---

	std::vector<std::string> varnames;
	
	// essential functions which must be defined by user
	virtual void set_size(std::array <double, dim> _x) = 0;
	virtual double init_density(std::array <double, dim> x, void * _env, double bf) = 0; 
	virtual std::array<double,dim> growthRate(std::array <double, dim> x, double t, void * _env) = 0;
	virtual double mortalityRate(std::array <double, dim> x, double t, void * _env) = 0;
	virtual double birthRate(std::array <double, dim> x, double t, void * _env) = 0;

	// optional functions which can be defined by user through overloads
	virtual ~IndividualBase(){
	}
	
	virtual void preCompute(std::array <double, dim> x, double t, void * _env){
	}

	virtual double establishmentProbability(double t, void  * _env){
		return 1;
	};

	// Functions related to additional state variables
	virtual void init_state(double t, void * _env){}

	virtual std::vector<double>::iterator set_state(std::vector<double>::iterator &it){
		return it;
	}
	virtual std::vector<double>::iterator get_state(std::vector<double>::iterator &it){
		return it;
	}
	virtual std::vector<double>::iterator get_rates(std::vector<double>::iterator &it){
		return it;
	} 

	// printing
	virtual void print(std::ostream &out = std::cout){}

	virtual void save(std::ostream& fout){}
	virtual void restore(std::istream& fin){}

};


#endif


