#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// Wrapper for GSL ODE Solver. Supplying Jacobian is not currently supported.

template <typename F> //, typename J>  
class ODE_Solver : public gsl_odeiv2_system {
	private:
	const F& _func;
//	const J& _jac;

	static int invoke_f(double t, const double y[], double dydt[], void *params) {
		return static_cast<ODE_Solver*>(params)->_func(t,y,dydt);
	}
//	static int invoke_j(double t, const double y[], double* dfdy, double dfdt[], void *params) {
//		return static_cast<ODE_Solver*>(params)->_jac(t,y,dfdy,dfdt);
//	}

	gsl_odeiv2_system * get(){
		return static_cast<gsl_odeiv2_system*>(this);
	}

	private:
	gsl_odeiv2_driver * d;
	
	public:
	double t;
	vector <double> y;

	public:
	ODE_Solver(const F& func, /*const J& jac,*/ size_t dim) : _func(func)/*, _jac(jac)*/ {
		function = &ODE_Solver::invoke_f;
		jacobian = NULL; //&ODE_Solver::invoke_j;
		dimension = dim;
		params=this;
		
		d = gsl_odeiv2_driver_alloc_y_new (get(), gsl_odeiv2_step_rkck, 1e-6, 1e-6, 0.0);
        t = 0;
	}

	
	~ODE_Solver(){
		gsl_odeiv2_driver_free(d);
	}
	
	inline void initialize(vector <double> &initial_state, double t0 = 0){
		y = initial_state;
		t = t0;
	}
	
	inline int step_to(double ti){
		int status = gsl_odeiv2_driver_apply (d, &t, ti, y.data());
		if (status != GSL_SUCCESS){
			cout << "Error in ODE_Solver, return value = (" << status << ")";
        }
		return status;
	}
};

//template <typename F>
//class ODE_Solver{
//	public:
//	gsl_odeiv2_driver * d;
//	ODE_Solver<F> *system;
//	
//	public:
//	vector <double> state;
//	double t;
//	
//	public:
//	ODE_Solver(const F func, size_t dim) {
//		system = new ODE_Solver<F>(func, dim);
//		d = gsl_odeiv2_driver_alloc_y_new (system->get(), gsl_odeiv2_step_rkck,1e-6, 1e-6, 0.0);
//	}
//	
//	~ODE_Solver(){
//		gsl_odeiv2_driver_free (d);
//	}
//	
//	inline void initialize(vector <double> &initial_state){
//		state = initial_state;
//		t = 0;
//	}
//	
//	inline int step_to(double ti){
////		cout << "Inside: " << state.size() << endl;
//		int status = gsl_odeiv2_driver_apply (d, &t, ti, state.data());
//		return status;
//	}
//}; 



//int
//func (double t, const double y[], double f[],
//      void *params)
//{
//  (void)(t); /* avoid unused parameter warning */
//  double mu = *(double *)params;
//  f[0] = y[1];
//  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
//  return GSL_SUCCESS;
//}

//int
//jac (double t, const double y[], double *dfdy,
//     double dfdt[], void *params)
//{
//  (void)(t); /* avoid unused parameter warning */
//  double mu = *(double *)params;
//  gsl_matrix_view dfdy_mat
//    = gsl_matrix_view_array (dfdy, 2, 2);
//  gsl_matrix * m = &dfdy_mat.matrix;
//  gsl_matrix_set (m, 0, 0, 0.0);
//  gsl_matrix_set (m, 0, 1, 1.0);
//  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
//  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
//  dfdt[0] = 0.0;
//  dfdt[1] = 0.0;
//  return GSL_SUCCESS;
//}

int main (void){
  double mu = 10;

  
  auto func_lambda = [&mu](double t, const double y[], double f[]) -> int {
	  f[0] = y[1];
	  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
	  return GSL_SUCCESS;
  };

  ODE_Solver <decltype(func_lambda)> ODESys(func_lambda, 2);
  
//  ODE_Solver <decltype(func_lambda)> odesys(func_lambda, 2);
  
  
//  gsl_odeiv2_system sys = {func, NULL, 2, &mu};

 
 
  double t = 0.0, t1 = 100.0;
//  double y[2] = { 1.0, 0.0 };

	vector <double> yvec = {1,0};
    ODESys.initialize(yvec, 0);


	ofstream fout("solution.txt");

  for (int i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
//      int status = gsl_odeiv2_driver_apply (ODESys.d, &t, ti, y);
		ODESys.step_to(ti);
//		cout << odesys.state.size() << endl;
      fout << ODESys.t << "\t" << ODESys.y[0] << "\t" << ODESys.y[1] << endl;
    }
	fout.close();

//  gsl_odeiv2_driver_free (d);
  return 0;
}



