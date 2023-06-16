#include "individual_base.h"
using namespace std;


IndividualBase::~IndividualBase(){
}

void IndividualBase::set_size(std::vector <double> _x){
}

void IndividualBase::preCompute(std::vector <double> x, double t, void * _env){
}

double IndividualBase::establishmentProbability(double t, void  * _env){
	return 1;
}

void IndividualBase::init_state(double t, void * _env){
}

vector<double>::iterator IndividualBase::set_state(vector<double>::iterator &it){
	return it;
}

vector<double>::iterator IndividualBase::get_state(vector<double>::iterator &it){
	return it;
}

vector<double>::iterator IndividualBase::get_rates(vector<double>::iterator &it){
	return it;
}

void IndividualBase::print(std::ostream& out){
}

void IndividualBase::save(std::ofstream &fout){
}

void IndividualBase::restore(std::ifstream &fin){
}
