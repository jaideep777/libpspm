
/// @param w            A function or function-object of the form `w(int i, double t)` that returns a double. It can be a lambda.
/// @param t            The current time (corresponding to the current physiological state), to be passed to `w`
/// @param xlow         The lower limit of the integral
/// @param species_id   The id of the species for which the integral should be computed
///
/// The integral is computed for different methods as follows:
///
/// \image html size_integral.png width=700cm 
///
/// `FMU:` \f$\quad I = \sum_{i=i_0}^J h_i w_i u_i\f$
/// 
/// `EBT:` \f$\quad I = \sum_{i=i_0}^J w_i N_i\f$, with \f$x_0 = x_b + \pi_0/N_0 \f$
/// 
/// `CM :` \f$\quad I = \sum_{i=i_0}^J h_i w_i (u_{i+1}+u_i)/2 \f$
/// 
/// where \f$i_0 = \text{argmax}(x_i \le x_{low})\f$, \f$h_i = x_{i+1}-x_i\f$, and \f$w_i = w(x_i)\f$.
/// 
/// Since the CM method uses the trapezoidal rule for integration, it includes the boundary cohort irrespective of the value of \f$x_0\f$.
 
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
	if (method == SOLVER_CM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		bool integration_completed = false;
		double I = 0;
		double u_hi = spp->getU(0); //(use_log_densities)? exp(spp->getU(0)) : spp->getU(0);
		double x_hi = spp->getX(0);
		double f_hi = w(0, t)*u_hi;
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		for (int i=1; i<spp->J; ++i){  // going from high ----> low
			double u_lo = spp->getU(i); 
			double x_lo = spp->getX(i);
			double f_lo = w(i, t)*u_lo;
	
			// This implementation allows stopping the integration exactly at xlow. TODO: Need to check if that actually makes sense
			//if (xlow < 0.01) cout << "x/w/u/f = " << x_lo << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_lo << "\n";
			if (x_lo <= xlow){
				// * interpolate to stop exactly at xlow
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
			
		if (spp->J == 1 || !integration_completed){  
			// boundary at xb
			double u0 = spp->get_boundary_u();
			double x_lo = spp->xb;
			double f_lo =  w(-1, t)*u0; // -1 is boundary cohort
			double f = (control.integral_interpolate)? f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo) : f_lo;  
			double x = (control.integral_interpolate)? xlow : x_lo;
			I += (x_hi-x)*(f_hi+f);
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
	
	else if (method == SOLVER_EBT){
		// set up cohorts to integrate
		// get (backup) pi0, N0 from last cohort 
		double   pi0  =  spp->getX(spp->J-1);
		double   N0   =  spp->getU(spp->J-1);

		// real-ize pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i){
		   	if (spp->getX(i) < xlow) break; // if X == xlow, we still include it in the intgral
			else I += w(i, t)*spp->getU(i);
		}
		
		// restore the original pi0-cohort
		spp->setX(spp->J-1, pi0);

		return I;
	}
	
	else{
		throw std::runtime_error("Unsupported solver method");
	}
}





//             _xm 
// Calculate _/ w(z,t)u(z,t)dz
//         xb
template<typename wFunc>
double Solver::integrate_x(wFunc w, double t, int species_id){
	Species_Base* spp = species_vec[species_id];

	if (method == SOLVER_FMU || method == SOLVER_IFMU){
		// integrate using midpoint quadrature rule
		double I=0;
		for (unsigned int i=0; i<spp->J; ++i){
			I += spp->h[i]*w(i, t)*spp->getU(i);  // TODO: Replace with std::transform after profiling
		}
		return I;
	}
	
	else if (method == SOLVER_EBT){
		// integrate using EBT rule (sum over cohorts)
		
		// backup pi0, N0 from last cohort 
		double   pi0  =  spp->getX(spp->J-1);
		double   N0   =  spp->getU(spp->J-1);

		// real-ize pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);

		// calculate integral
		double I = 0;
		for (int i=0; i<spp->J; ++i) I += w(i, t)*spp->getU(i);
		
		// restore the original pi0-cohort
		spp->setX(spp->J-1, pi0);

		return I;
	}
	
	if (method == SOLVER_CM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the beginning, so x will be ascending
		double I = 0;
		double u_hi = spp->getU(0); //(use_log_densities)? exp(spp->getU(0)) : spp->getU(0);
		double x_hi = spp->getX(0);
		double f_hi = w(0, t)*u_hi;
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";

		for (int i=1; i<spp->J; ++i){
			double u_lo = spp->getU(i); //(use_log_densities)? exp(spp->getU(i)) : spp->getU(i);
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
	
	else{
		throw std::runtime_error("Unsupported solver method");
	}
}


