#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>

template<class Ind>
class Cohort : public Ind {
	public:
	double x = 1;
	double u = 2;
	int id = 0;	

	void print(){
		std::cout << "id = " << id << ", x = " << x << ", u = " << u << "\n"; 
	}	
};

#endif
