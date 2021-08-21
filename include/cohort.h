#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>

template<class Ind>
class Cohort : public Ind {
	public:
	double x = -999;
	double u = -999;
	int id = 0;	

	double birth_time;
	bool remove = false;

	void print_xu(){
		std::cout << id << "\t" << x << "\t" << u << "\t"; 
	}

	void set_size(double _x){
		x = _x;
		Ind::set_size(x);
	}	
};

#endif
