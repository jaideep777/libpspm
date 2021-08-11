#include <iostream>
#include <cohort.h>

class Individual {
	public:
	double lma = 30;
};

int main(){
	Cohort<Individual> C;

	std::cout << " " << C.lma << " " << C.x << " " << C.u << " " << C.id << "\n";
}
