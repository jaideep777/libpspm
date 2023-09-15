#include <iostream>
#include <vector>
#include <environment_base.h>
#include <individual_base.h>
#include <cohort.h>
#include <string_view>
using namespace std;


class LightEnv : public EnvironmentBase{
	public:
	double E = 0.95;

	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
		E = 0.95 + 0.05*t;
	}

};



class Plant : public IndividualBase<1>{
	public:
	double lma = 30;

	double height = -99;
	double crown_area = -99;
	double root_mass = -99;

	vector<string> varnames = {"lma|", "ht", "cr", "root"};

	Plant(){
		lma = 10;
	}

	Plant(double _lma){
		lma = _lma;
	}
	
	void set_size(array<double,1> x){
		height = x[0];
		crown_area = height*height;
	}

	void preCompute(array<double,1> x, double t, void * env){
	}

	array<double,1> growthRate(array<double,1> x, double t, void * env){
		cout << "in g: " << x[0] << " " << t << " " << ((LightEnv*)env)->E << "\n";
		return {x[0]*((LightEnv*)env)->E};
	}
	double mortalityRate(array<double,1> x, double t, void * env){
		return -0.5;
	}
	double birthRate(array<double,1> x, double t, void * env){
		return 1;
	}
	
	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	double init_density(array<double,1> x, void * _env, double bf){
		return 5/(x[0]+0.5);
	}

	void init_state(double t, void * env){
		root_mass = 0;
	}

	vector<double>::iterator set_state(vector<double>::iterator &it){
		root_mass = *it++;
		return it;
	}

	vector<double>::iterator get_state(vector<double>::iterator &it){
		*it++ = root_mass;
		return it;
	}

	vector<double>::iterator get_rates(vector<double>::iterator &it){
		*it++ = -0.1*(root_mass-1);
		return it;
	}

	void print(std::ostream& out = std::cout){
		out << std::setw(6) << std::setprecision(4) 
		    << lma << "\t" << height << "\t" << crown_area << "\t" << root_mass << "\t";
	}
	
	void save(std::ofstream& fout){
	}
	void restore(std::ifstream& fin){
	}
	
};


class Insect : public IndividualBase<2>{
	public:
	double w; // weight
	double e; // energy
	double wingspan = 10;
	
	double m = -1, f = -1;
	array<double,2> g;

	vector<string> varnames = {"wing", "f"};
	
	double init_density(array<double,2> x, void * _env, double bf){
		return x[0]*x[1]/10;
	}

	void set_size(array<double,2> x){
		w = x[0];
		e = x[1];
		wingspan = 2*w;	
	}

	void preCompute(array<double,2> x, double t, void * env){
		LightEnv * le = static_cast<LightEnv*>(env);
		f = 0.1*w;
		g = {w*0.1, le->E*0.1};
		m = w/e;
	}

	array<double,2> growthRate(array<double,2> x, double t, void * env){
		return g;
	}
	double mortalityRate(array<double,2> x, double t, void * env){
		return m;
	}
	double birthRate(array<double,2> x, double t, void * env){
		return f;
	}

	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	void print(std::ostream& out = std::cout){
		out << std::setw(10) << std::setprecision(4) 
		    << w << '\t' << e << '\t'
			<< wingspan << '\t' << f << "\t";
	}
};


template <class Model>
class TestSpecies{
	public:
	Cohort<Model> c;
	auto getX(int i){
		return c.x;
	}
};


template <typename T>
constexpr auto type_name() {
  std::string_view name, prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void)";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

int main(){

	LightEnv E;

	Plant P;
	P.set_size({2});
	cout << "P: "; P.print(); cout << "\n";

	P.lma = 70;
	Cohort<Plant> C1;
	C1 = Cohort<Plant> (P);
	cout << "C1: "; C1.print(); cout << "\n";

	Insect I1;
	I1.set_size({10, 90});
	cout << "I: "; I1.print(); cout << "\n";

	Cohort<Insect> C2(I1);
	cout << "C2: "; C2.print(); cout << "\n";

	C2.set_size(array<double,2>{40,60});
	cout << "C2: "; C2.print(); cout << "\n";

	E.computeEnv(1, nullptr, vector<double>().begin(), vector<double>().begin());
	cout << "Env @ t = 1: " << E.E << '\n';
	C2.preCompute(array<double,2>{0,0}, 1, &E);
	cout << "Insect g/m/f: " << C2.g << " / " << C2.m << " / " << C2.f << '\n';

	C2.save(cout, 0);

	TestSpecies<Insect> Si;
	Si.c = C2;
	cout << "Cohort inside Insect species has state: " << Si.getX(0) << '\n';
	cout << "Cohort inside Insect species has state type: " << type_name<decltype(Si.getX(0))>() << '\n';

	TestSpecies<Plant> Sp;
	Sp.c = C1;
	cout << "Cohort inside Plant species has state: " << Sp.getX(0) << '\n';
	cout << "Cohort inside Plant species has state type: " << type_name<decltype(Sp.getX(0))>() << '\n';




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



