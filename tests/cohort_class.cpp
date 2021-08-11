#include <iostream>
#include <cohort.h>
#include <species.h>
#include <vector>
using namespace std;


class Individual {
	public:
	double lma = 30;

	double height;
	double crown_area;
	double root_mass;

	int nstate = 1;

	double set_height(double x){
		height = x;
		crown_area = x*x;
	}

	vector<double>::iterator set_state(vector<double>::iterator &it){
		root_mass = *it++;
		return it;
	}

	vector<double>::iterator get_rates(vector<double>::iterator &it){
		*it++ = -0.1*(root_mass-1);
		return it;
	}
};

int main(){
	Cohort<Individual> C;

	std::cout << " " << C.lma << " " << C.x << " " << C.u << " " << C.id << "\n";
}



