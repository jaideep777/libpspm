// -*-c++-*-
#ifndef PLANT_PLANT_ODE_SOLVER_H_
#define PLANT_PLANT_ODE_SOLVER_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <vector>
#include <chrono>


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
	double abs_tol = 1e-4;
	double rel_tol = 1e-4;
	double initial_step_size = 1e-6;
	
	public:
	double t;
	std::vector <double> y;

	public:
	ODE_Solver(const F& func, /*const J& jac,*/ size_t dim) : _func(func)/*, _jac(jac)*/ {
		function = &ODE_Solver::invoke_f;
		jacobian = NULL; //&ODE_Solver::invoke_j;
		dimension = dim;
		params=this;
		
		d = gsl_odeiv2_driver_alloc_y_new (get(), gsl_odeiv2_step_rkck, initial_step_size, abs_tol, rel_tol);
//		gsl_odeiv2_driver_set_hmin(d, 1e-6);
        t = 0;
	}

	
	~ODE_Solver(){
		gsl_odeiv2_driver_free(d);
	}
	

	inline void initialize(std::vector <double> &initial_state, double t0 = 0){
		y = initial_state;
		t = t0;
	}
	

	inline int step_to(double ti){
//		auto t1 = std::chrono::steady_clock::now();
		int status = gsl_odeiv2_driver_apply (d, &t, ti, y.data());
//		auto t2 = std::chrono::steady_clock::now();
		  
		if (status != GSL_SUCCESS){
			std::cerr << "Error in ODE_Solver, return value = (" << status << ")";
        }
//		std::cout << "Stepped with " << d->n << " steps, step size = " << d->h << " [" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/1e3 << " usec]" << std::endl;
		return status;
	}
	

	inline void resize(size_t newdim){
		dimension = newdim;
		gsl_odeiv2_driver_free(d);
		d = gsl_odeiv2_driver_alloc_y_new (get(), gsl_odeiv2_step_rkck, initial_step_size, abs_tol, rel_tol);
//		gsl_odeiv2_driver_set_hmin(d, 1e-6);
	}


//	inline int step_to_euler(double ti){
//		double dt = ti-t;
//		std::vector <double> f(y.size());
//		_func(t, y.data(), f.data());
//		for (int i=0; i<f.size(); ++i) std::cout << f[i] << " "; 
//		std::cout << std::endl;
//		for (size_t i=0; i<y.size(); ++i) y[i] += f[i]*dt;
//		return 0;
//		t = ti;
//	}
	
};



#endif

