#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>

template<class Ind>
class Cohort : public Ind {
	public:
	double x = -999;
	double u = -999;
	int id = 0;	

	double birth_time = 0;
	bool remove = false;

	// Default constructor simply calls Individual's default ctor
	// NOTE: This means user's Ind class will need ctor with no arguments
	Cohort() : Ind() {
	}

	// Construct a cohort from Individual using copy constructor of Individual
	Cohort(const Ind& _ind) : Ind(_ind){
	}

	void print_xu(std::ostream &out = std::cout){
		out << birth_time << "\t" << x << "\t" << u << "\t"; 
	}

	void print(std::ostream &out = std::cout){
		print_xu(out);
		Ind::print(out);
	}

	void set_size(double _x){
		x = _x;
		Ind::set_size(x);
	}
};

#endif

