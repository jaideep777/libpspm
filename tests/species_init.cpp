#include <iostream>
#include <vector>
#include <cohort.h>
#include <species.h>
#include <cmath>
using namespace std;

#include "test_model_plant_insect.h"

template<class Model>
int check_species_states(Species<Model>& spp, const vector<double>& expected){
	for (int i=0, k=0; i<spp.xsize(); ++i){
		auto& c = spp.getCohort(i);
		for (auto& x : c.x){
			if (fabs(x - expected[k]) > 1e-6) return 1;
			++k; 
		}
	}
	return 0;
}

template<class T>
int check(const vector<T>& v1, const vector<T>& v2){
	if (v1.size() != v2.size()) return 1;

	for (int i=0; i<v1.size(); ++i){
		if (fabs(v1[i]-v2[i]) > 1e-6) return 1;
	}
	return 0;
}

int check(double v1, double v2){
	if (fabs(v1-v2) > 1e-6) return 1;
	else return 0;
}


int main(){

	int nerrors = 0;

	LightEnv E;

	Plant P;
	P.set_size({25});
	cout << "P: "; P.print(); cout << "\n";

	P.lma = 70;
	Cohort<Plant> C1;
	C1 = Cohort<Plant> (P);
	C1.set_size({30});
	cout << "C1: "; C1.print(); cout << "\n";

	Insect I1;
	I1.set_size({10, 90});
	cout << "I: "; I1.print(); cout << "\n";

	Cohort<Insect> C2(I1);
	cout << "C2: "; C2.print(); cout << "\n";

	C2.set_size({40,60});
	cout << "C2: "; C2.print(); cout << "\n";

	E.computeEnv(1, nullptr, vector<double>().begin(), vector<double>().begin());
	cout << "Env @ t = 1: " << E.E << '\n';
	C2.preCompute(1, &E);
	cout << "Insect g/m/f: " << C2.g << " / " << C2.m << " / " << C2.f << '\n';
	C2.print(); cout << '\n';

	C2.save(cout, 0);

	Species<Insect> Si(I1); // This constructor does not necessarily set x
	Si.set_xb({0.5, 0.01});
	Si.print();

	Cohort<Insect> C3; C3.set_size({20, 0.01}); C3.u = 1.1;
	Cohort<Insect> C4; C4.set_size({20, 10});   C4.u = 1.3;
	Cohort<Insect> C5; C5.set_size({0.5, 5});   C5.u = 2.5;
	Si.addCohort(C3);
	Si.addCohort(C4);
	Si.addCohort(C5);
	Si.print();

	Species<Plant> Sp(P);
	Sp.set_xb({0.35});
	Sp.print();	

	Cohort<Plant> Cp3; Cp3.set_size({30}); Cp3.u = 5.1;
	Cohort<Plant> Cp4; Cp4.set_size({31}); Cp4.u = 3.1;
	Cohort<Plant> Cp5; Cp5.set_size({32}); Cp5.u = 1.1;
	Cohort<Plant> Cp6; Cp6.set_size({32.1}); Cp6.u = 0.1;

	Sp.addCohort(Cp3);
	Sp.addCohort(Cp4);
	Sp.addCohort(Cp5);
	Sp.addCohort(Cp6);

	Sp.print();
	nerrors += check_species_states(Sp, {30, 31, 32, 32.1});

	cout << "Sort cohorts by body size\n-----------------------\n";
	Si.sortCohortsDescending(0);
	Si.print();
	nerrors += check_species_states(Si, {20, 0.01, 20, 10, 0.5, 5});

	cout << "Sort cohorts by energy reserve\n-----------------------\n";
	Si.sortCohortsDescending(1);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 20, 0.01});

	// Now insert a EBT boundary cohort and check that it can be skipped in sorting
	Cohort<Insect> Cb; Cb.set_size({100, 0.02}); Cb.u = 0;
	Si.addCohort(Cb);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 0.5, 5, 20, 0.01, 100, 0.02});

	cout << "Sort cohorts by body size, exclude BC\n-----------------------\n";
	Si.sortCohortsDescending(0, 1);
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 20, 0.01, 0.5, 5, 100, 0.02});


	cout << "Set state of a specific cohort\n-----------------------\n";
	Si.setX(1, {50, 20});
	Si.print();
	nerrors += check_species_states(Si, {20, 10, 50, 20, 0.5, 5, 100, 0.02});
	
	auto g1 = Si.growthRate(2, 1, &E); // cohort 2 {0.5, 5}
	cout << g1 << '\n';
	nerrors += check(g1, {0.05, 0.1});
	nerrors += check(Si.mortalityRate(2, 1, &E), 0.5/5);
	nerrors += check(Si.birthRate(2, 1, &E), 0.05);

	nerrors += check(Sp.growthRate(2, 1, &E), {32});

	return nerrors;

// 	Species<Plant> Sv(array<double,1> {1,2,3,4,5});
// 	Sv.setX(0, 1.5);
// 	Sv.setX(2, 3.5);
// 	Sv.print();

// 	Species<Plant> S(P);  // create species S with prototype P

// 	Species<Insect> I(array<double,1> {1.1,2.1,3.1});

// 	Species_Base * S1 = &S;
// 	//S1->print();
// 	cout << "S size: " << S1->get_maxSize() << "\n";

// 	Species_Base * S2 = &I;
// 	//S2->print();
// 	cout << "I size: " << S2->get_maxSize() << "\n";

// 	LightEnv env;

// 	Solver sol(SOLVER_EBT);
// 	sol.setEnvironment(&env);
// 	sol.addSpecies(array<double,1> {1,2,3,4,5}, &S, 1, 1);
// 	sol.addSpecies(array<double,1> {1.1,2.1,3.1}, &I, 0, 1);
// 	sol.resetState();
// 	sol.initialize();
// 	//sol.addCohort_EBT();
// 	sol.print();	

// 	I.setX(2, 0.05);
// 	I.setU(2, 0.2);
// 	sol.print();	
// 	cout << "Insect BF = " << sol.calcSpeciesBirthFlux(1,0) << "\n";
	
// //	sol.preComputeSpecies(1,0);
// 	sol.print();	
// 	cout << "Insect BF (after precompute) = " << sol.calcSpeciesBirthFlux(1,0) << "\n";

// 	Cohort<Plant> C;
// 	C.print(); cout << "\n";
// 	C.set_size(10);
// 	C.print(); cout << "\n";
// 	cout << C.growthRate(C.height, 0, sol.env) << "\n";

}



