
/// @param w            A function or function-object of the form `w(int i, double t)` that returns a double. 
///                     It can be a lambda. This function should access the 'i'th cohort and compute the weight 
///                     from the cohort's properties. The function `w` should be able access to the `Solver` in
///                     order to access cohorts. 
/// @param t            The current time (corresponding to the current physiological state). This will be passed 
///                     to `w` to allow the weights to be a direct function of time.
/// @param xlow         The lower limit of the integral
/// @param species_id   The id of the species for which the integral should be computed
///
/// \image html size_integral.png width=700cm 
///
/// the computation of this integral depends on the solver method. For different solvers, the integral is defined as follows:
/// 
/// `FMU:` \f$\quad I = \sum_{i=i_0}^J h_i w_i u_i\f$
/// 
/// `EBT:` \f$\quad I = \sum_{i=i_0}^J w_i N_i\f$, with \f$x_0 = x_b + \pi_0/N_0\f$
/// 
/// `CM :` \f$\quad I = \sum_{i=i_0}^J h_i  (w_{i+1}u_{i+1}+w_i u_i)/2\f$
/// 
/// where \f$i_0 = \text{argmax}(x_i \le x_{low})\f$, \f$h_i = x_{i+1}-x_i\f$, and \f$w_i = w(x_i)\f$. 
/// 
/// If interpolation is turned on, \f$h_{i_0}=x_{i_0+1}-x_{low}\f$, whereas \f$u(x_{low})\f$ is set 
/// to \f$u(x_{i_0})\f$ in FMU and calculated by bilinear interpolation in CM (See Figure). 
/// Interpolation does not play a role in EBT.
/// 
/// In the CM method, the density of the boundary cohort is obtained from the boundary condition 
/// of the PDE: \f$u_b=B/g(x_b)\f$, where \f$B\f$ is the input flux of newborns. In the current implementation 
/// of the CM method, \f$B\f$ must be set to a constant. In future implementations which may allow 
/// real-time calculation of \f$B\f$, \f$u_b\f$ must be calculated recurvisely by solving \f$u_b = B(u_b)/g(x_b)\f$, 
/// where \f$B(u_b) = \int_{x_b}^{x_m} f(x)u(x)dx\f$.

//             _xm 
// Calculate _/ w(z,t)u(z)dz
//         xlow
// implementation from orig plant model	
// ----
// I += (x_hi - x_lo) * (f_hi + f_lo);
// x_hi = x_lo;
// f_hi = f_lo;
// if (x_lo < xlow) break;
// ---- 
template<typename wFunc>
double Solver::integrate_wudx_above(wFunc w, double t, double xlow, int species_id){

	Species_Base* spp = species_vec[species_id];

	// cohorts are inserted at the end, hence sorted in descending order
	// FIXME: should there be an optional "sort_cohorts" control parameter? Maybe some freak models are such that cohorts dont remain in sorted order?
	if (method == SOLVER_CM || method == SOLVER_ICM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		bool integration_completed = false;
		double I = 0;
		double u_hi = spp->getU(0); 
		double x_hi = spp->getX(0);
		double f_hi = w(0, t)*u_hi;
		
		for (int i=1; i<spp->J; ++i){  // going from high ----> low
			double u_lo = spp->getU(i); 
			double x_lo = spp->getX(i);
			double f_lo = w(i, t)*u_lo;
	
			if (x_lo <= xlow){  // check if xlow falls within the current interval. 
				// * if interpolation is on, stop exactly at xlow, else take the full interval
				double f = (control.integral_interpolate)? f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo) : f_lo;  
				double x = (control.integral_interpolate)? xlow : x_lo;
				I += (x_hi - x) * (f_hi + f);
				integration_completed = true;
			}
			else{
				I += (x_hi - x_lo) * (f_hi + f_lo);
			}
		
			x_hi = x_lo;
			f_hi = f_lo;
			if (integration_completed) break;
		}
		
		//if (integration_completed) assert(x_hi < xlow);
		//if (!integration_completed) assert(x_hi >= xlow || spp.J == 1);
		
		// if integration has not completed, continue to the boundary interval	
		if (spp->J == 1 || !integration_completed){  
			// boundary at xb
			double u0 = spp->get_boundary_u();
			double x_lo = spp->xb;
			if (x_hi > x_lo){ 
				double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
				double f = (control.integral_interpolate)? f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo) : f_lo;  
				double x = (control.integral_interpolate)? xlow : x_lo;
				I += (x_hi-x)*(f_hi+f);
			}
		}
		
		return I*0.5;
	}

	else if (method == SOLVER_FMU || method == SOLVER_IFMU){
		//if (xlow < spp->xb) throw std::runtime_error("integral lower bound must be >= xb");
		if (xlow > spp->x[spp->J]) return 0;  // if xlow is above the maximum size, integral is 0

		// integrate using midpoint quadrature rule
		double I=0;
		for (int i=spp->J-1; i>=0; --i){  // in FMU, cohorts are sorted ascending
			if (spp->x[i] <= xlow){ // check if last interval is reached 
				double h = (control.integral_interpolate)?  spp->x[i+1]-xlow  :  spp->h[i];
				I += h * w(i,t) * spp->getU(i);
				break;
			}
			else{
				I += spp->h[i] * w(i,t) * spp->getU(i);
			}
		}
		
		return I;
	}
	
	else if (method == SOLVER_EBT || method == SOLVER_IEBT){
		// set up cohorts to integrate
		realizeEbtBoundaryCohort(spp);

		// sort cohorts, but skip pi0-cohort
		spp->sortCohortsDescending(1); // FIXME: Add a label to EBT boundary cohort and assert that it is always at J-1

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){  // in EBT, cohorts are sorted descending
		   	if (spp->getX(i) < xlow) break; // if X == xlow, we still include it in the intgral
			else I += w(i, t)*spp->getU(i);
		}
		
		// restore the original pi0-cohort
		restoreEbtBoundaryCohort(spp);

		return I;
	}
	
	else if (method == SOLVER_ABM){
		// sort cohorts descending - this is important in ABM because traits may vary
		spp->sortCohortsDescending(); 
		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){       // in ABM, cohorts are sorted descending
		   	if (spp->getX(i) < xlow) break; // if X == xlow, we still include it in the intgral
			else I += w(i, t)*spp->getU(i);
		}
		
		return I;
	}
	
	else{
		throw std::runtime_error("Unsupported solver method");
	}
}





//             _xm 
// Calculate _/ w(z,t)u(z,t)dz
//         xb
// This will eventually be replaced with a call to integrate_wudx_above
template<typename wFunc>
double Solver::integrate_x(wFunc w, double t, int species_id){
	Species_Base* spp = species_vec[species_id];

	if (method == SOLVER_FMU || method == SOLVER_IFMU){
		// integrate using midpoint quadrature rule
		double I=0;
		for (unsigned int i=0; i<spp->J; ++i){
			I += spp->h[i]*w(i, t)*spp->getU(i);  // TODO: Replace with std::transform after profiling
			//std::cout << "integral: " << w(i, t) << " * " << spp->getU(i) << std::endl;
		}
		return I;
	}
	
	else if (method == SOLVER_EBT || method == SOLVER_IEBT){
		// integrate using EBT rule (sum over cohorts)
		realizeEbtBoundaryCohort(spp);

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i) I += w(i, t)*spp->getU(i);
		
		// restore the original pi0-cohort
		restoreEbtBoundaryCohort(spp);

		return I;
	}
	
	if (method == SOLVER_CM || method == SOLVER_ICM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		double I = 0;
		double u_hi = spp->getU(0); 
		double x_hi = spp->getX(0);
		double f_hi = w(0, t)*u_hi;

		for (int i=1; i<spp->J; ++i){
			double u_lo = spp->getU(i); 
			double x_lo = spp->getX(i);
			double f_lo = w(i, t)*u_lo;
	
			I += (x_hi - x_lo) * (f_hi + f_lo);
			x_hi = x_lo;
			f_hi = f_lo;
		}
		
		// boundary at xb
		double u0 = spp->get_boundary_u();
		double x_lo = spp->xb;
		double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
		I += (x_hi-x_lo)*(f_hi+f_lo);
		
		return I*0.5;
	}
	
	else if (method == SOLVER_ABM){
		// calculate integral. Sorting is not required here because all cohorts will be touched anyway
		double I = 0;
		for (int i=0; i<spp->J; ++i){       // in ABM, cohorts are sorted descending
			I += w(i, t)*spp->getU(i);
		}
		
		return I;
	}
	
	else{
		throw std::runtime_error("Unsupported solver method");
	}
}


