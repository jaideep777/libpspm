template<typename wFunc>
double Solver::integrate_wudx_above(wFunc w, double t, double xlow, int species_id){

	Species_Base* spp = species_vec[species_id];

	// cohorts are inserted at the end, hence sorted in descending order - FIXME: should there be an optional "sort_cohorts" control parameter? Maybe some freak models are such that cohorts dont remain in sorted order?
	if (method == SOLVER_CM){
		//std::cout << "J = " << spp.J << ", dist = " << iset.dist << std::endl; 
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the end, so x will be descending
		bool integration_completed = false;
		double I = 0;
		double u_hi = spp->getU(0); //(use_log_densities)? exp(spp->getU(0)) : spp->getU(0);
		double x_hi = spp->getX(0);
		double f_hi = w(0, x_hi, t)*u_hi;
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		//--itx; --itu;
		for (int i=1; i<spp->J; ++i){
			double u_lo = spp->getU(i); //(use_log_densities)? exp(spp->getU(i)) : spp->getU(i);
			double x_lo = spp->getX(i);
			double f_lo = w(i, x_lo, t)*u_lo;
	
			//// implementation from orig plant model	
			//I += (x_hi - x_lo) * (f_hi + f_lo);
			//x_hi = x_lo;
			//f_hi = f_lo;
			//if (x_lo < xlow) break;
			// ---- 

			// This implementation allows stopping the integration exactly at xlow. TODO: Need to check if that actually makes sense
			//if (xlow < 0.01) cout << "x/w/u/f = " << x_lo << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_lo << "\n";
			if (x_lo < xlow){
				// * interpolate to stop exactly at xlow
				//double f = f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo);  
				//double x = xlow;
				// * include integration up to next cohort
				double f = f_lo;
				double x = x_lo;

				I += (x_hi - x) * (f_hi + f);
				x_hi = x_lo;  // - these two probably not needed
				f_hi = f_lo;  // - 
				integration_completed = true;
				break;
			}
			else{
				I += (x_hi - x_lo) * (f_hi + f_lo);
				x_hi = x_lo;
				f_hi = f_lo;
			}
		}
		
		//if (integration_completed) assert(x_hi < xlow);
		//if (!integration_completed) assert(x_hi >= xlow || spp.J == 1);
			
		// if (spp.J == 1 || (f_hi > 0)){  // <-- this is from original plant model
		// .. (f_hi > 0) condition works for plant model but may not work generically. 
		// Instead use an explicit flag to mark completion 
		//if (spp.J == 1 || !integration_completed){  
			////double g = spp.mod->growthRate(spp.xb, t, env);
			////double u0 = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(t, env)/g  :  0; 
			//// if (x_b < xlow){ // interpolate. FIXME: Need the interpolation condition here too. 
			//double u0 = spp.u0_save;
			//double x_lo = spp.xb;	
			//double f_lo =  w(x_lo, t)*u0;
			//I += (x_hi-x_lo)*(f_hi+f_lo);
		//}
		
		return I*0.5;
	}

	//else if (method == SOLVER_FMU){
		//// integrate using midpoint quadrature rule
		//double I=0;
		//auto iset = spp.get_iterators(S);
		//auto &itu = iset.get("u");
		//iset.begin();
		//double * U = &(*itu);

		//for (int i=spp.J-1; i>=0; --i){
			////cout << "Enter: " << i << endl;
			////I += spp.h[i]*w(spp.X[i], t)*U[i];  // TODO: Replace with std::transform after profiling
			//double f = w(spp.X[i],t)*U[i];
			////std::cout << "f = " << f << " " << spp.x[i] << " " << xlow << " " << spp.h[i] << std::endl;
			//if (spp.x[i] < xlow){
				////I += (spp.x[i+1]-xlow) * f; // interpolating the last interval
				//I += spp.h[i] * f; // including full last interval
				//break;
			//}
			//else{
				//I += spp.h[i] * f;
				////std::cout << "Here: " << i << std::endl;
			//}
		//}
		////std::cout << "Here" << std::endl;
		//return I;
	//}
	
	//else if (method == SOLVER_EBT){
		//// integrate using midpoint quadrature rule
		//double I=0;
		//auto iset = spp.get_iterators(S);
		//auto &itu = iset.get("u");
		//auto &itx = iset.get("X");
		
		//double pi0 = *itx;
		//double N0  = *itu;

		//iset.rbegin();

		//for (int i=spp.J-1; i>=1; --i, --iset){  // iterate over cohorts except boundary cohort
			//if (*itx < xlow) break;
			
			//double f = w(*itx,t) * (*itu);
			////std::cout << "f = " << f << " " << spp.x[i] << " " << xlow << " " << spp.h[i] << std::endl;
			//I += f;
		//}
		
		//double x0 = spp.xb + pi0/(N0+1e-12);  
		//if (xlow < x0) I += w(x0, t)*N0;	 
		
		////std::cout << "Here" << std::endl;
		//return I;
	//}
	
	//else{
		//std::cout << "Only CM is implemented\n";
		//return 0;
	//}
}






template<typename wFunc>
double Solver::integrate_x(wFunc w, double t, int species_id){
	Species_Base* spp = species_vec[species_id];

	//if (method == SOLVER_FMU){
		//// integrate using midpoint quadrature rule
		//double I=0;
		//double * U = &(*itu);
		//for (unsigned int i=0; i<spp.J; ++i){
			//I += spp.h[i]*w(spp.X[i], t)*U[i];  // TODO: Replace with std::transform after profiling
		//}
		//return I;
	//}
	
	//else if (method == SOLVER_EBT){
		//// integrate using EBT rule (sum over cohorts)
		//iset.begin();
		//double   pi0  =  *itx;
		//double   N0   =  *itu;
		//++iset; // skip boundary cohort

		//double x0 = spp.xb + pi0/(N0+1e-12); 
		//double I = w(x0, t)*N0;
		//for (; !iset.end(); ++iset) I += w(*itx, t)*(*itu);
		
		//return I;
	//}
	
	if (method == SOLVER_CM){
		// integrate using trapezoidal rule. Below modified to avoid double computation of w(x)
		//double I = 0;
		//for (iset.begin(); iset.dist < iset.size-1; ++iset){
		//    double unext = (use_log_densities)? exp(*(itu+1)) : *(itu+1);
		//    double unow  = (use_log_densities)? exp(*itu) : *itu;
		//    I += (*(itx+1)-*itx)*(w(*(itx+1), t)*unext + w(*itx, t)*unow);
		//}
		//return I*0.5;
		
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the beginning, so x will be ascending
		double I = 0;
		double u_hi = spp->getU(0); //(use_log_densities)? exp(spp->getU(0)) : spp->getU(0);
		double x_hi = spp->getX(0);
		double f_hi = w(0, x_hi, t)*u_hi;
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		//--itx; --itu;
		for (int i=1; i<spp->J; ++i){
			double u_lo = spp->getU(i); //(use_log_densities)? exp(spp->getU(i)) : spp->getU(i);
			double x_lo = spp->getX(i);
			double f_lo = w(i, x_lo, t)*u_lo;
	
			I += (x_hi - x_lo) * (f_hi + f_lo);
			x_hi = x_lo;
			f_hi = f_lo;
		}
		
		// boundary at xb
		double u0 = spp->get_boundary_u();
		double x_lo = spp->xb;
		double f_lo =  w(spp->J, x_lo, t)*u0; // FIXME: spp->J is incorrect, how to specify boundary cohort?
		I += (x_hi-x_lo)*(f_hi+f_lo);
		
		return I*0.5;
	}
	
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


