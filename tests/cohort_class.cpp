#include <iostream>
#include <cohort.h>
#include <species.h>
#include <vector>
using namespace std;


class Plant {
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


class Insect {
	public: 
};


int main(){
	Cohort<Plant> C;
	
	Species<Plant> S(vector<double> {1,2,3,4,5});
	S.print();

	Species<Insect> I(vector<double> {1.1,2.1,3.1});

	Species_Base * S1 = &S;
	S1->print();
	cout << "S size: " << S1->get_maxSize() << "\n";

	Species_Base * S2 = &I;
	S2->print();
	cout << "I size: " << S2->get_maxSize() << "\n";

	std::cout << " " << C.lma << " " << C.x << " " << C.u << " " << C.id << "\n";
}



