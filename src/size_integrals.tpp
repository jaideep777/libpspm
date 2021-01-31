template<class Model, class Environment>
template<typename wFunc>
double Solver<Model, Environment>::integrate_wudx_above(wFunc w, double t, double xlow, vector<double>&S, int species_id){
	//cout << " | " <<  t << " " << mod->evalEnv(0,t) << " ";
	//if (method == SOLVER_FMU){
	//}
	//else if (method == SOLVER_EBT){
	//}
	//cout << "Begin integrate: xsize = " << xsize() << "(" << S[0] << ", " << S[xsize()-1] << "), xlow = " << xlow << endl;
	Species<Model> &spp = species_vec[species_id];
	auto iset = spp.get_iterators(S);
	auto &itx = iset.get("X");
	auto &itu = iset.get("u");
	iset.rbegin();
	//std::cout << "J = " << spp.J << ", dist = " << iset.dist << std::endl; 

	if (method == SOLVER_CM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the beginning, so x will be ascending
		bool converged = false;
		double I = 0;
		double u = (use_log_densities)? exp(*itu) : *itu;
		double x_hi = *itx;
		double f_hi = w(*itx, t)*u;
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		--iset; //--itx; --itu;
		for (int i=0; i<spp.J-1; ++i){
			double u = (use_log_densities)? exp(*itu) : *itu;
			double x_lo = *itx;
			double f_lo = w(*itx,t)*u;
			//if (xlow < 0.01) cout << "x/w/u/f = " << x_lo << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_lo << "\n";
			--iset; //--itx; --itu;
			if (x_lo < xlow){
				double f = f_lo; //f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo);  // FIXME: these should stop at the interpolating point
				double x = x_lo; //xlow;
				I += (x_hi-x) * (f_hi + f);
				converged = true;
				break;
			}
			else{
				I += (x_hi - x_lo) * (f_hi + f_lo);
			}
			x_hi = x_lo;
			f_hi = f_lo;
		}
		
		//if (spp.J == 1 || (f_hi > 0 && !converged)){  // f_hi condition causes interpolator to break. Need to check
		//    double x_lo = spp.xb;	// FIXME: should be max(spp.xb, xlow)
		//    double g = spp.mod->growthRate(spp.xb, t, env);
		//    double u0 = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(t, env)/g  :  0; 
		//    double f_lo =  w(spp.xb, t)*u0;
		//    //double f_lo = f_hi*2;
		//    I += (x_hi-x_lo)*(f_hi+f_lo);
		//}
		return I*0.5;
	}
	else{
		std::cout << "Only CM is implemented\n";
		return 0;
	}
}






template<class Model, class Environment>
template<typename wFunc>
double Solver<Model,Environment>::integrate_x(wFunc w, double t, vector<double>&S, int species_id){
	Species<Model> &spp = species_vec[species_id];
	auto iset = spp.get_iterators(S);
	auto &itx = iset.get("X");
	auto &itu = iset.get("u");

	//cout << " | " <<  t << " " << mod->evalEnv(0,t) << " ";
	if (method == SOLVER_FMU){
		// integrate using midpoint quadrature rule
		double I=0;
		double * U = &(*itu);
		for (unsigned int i=0; i<spp.J; ++i){
			I += spp.h[i]*w(spp.X[i], t)*U[i];  // TODO: Replace with std::transform after profiling
		}
		return I;
	}
	else if (method == SOLVER_EBT){
		// integrate using EBT rule (sum over cohorts)
		iset.begin();
		double   pi0  =  *itx;
		double   N0   =  *itu;
		++iset;

		double x0 = spp.xb + pi0/(N0+1e-12); 
		double I = w(x0, t)*N0;
		for (; !iset.end(); ++iset) I += w(*itx, t)*(*itu);
		
		return I;
	}
	else if (method == SOLVER_CM){
		// integrate using trapezoidal rule TODO: Modify to avoid double computation of w(x)
		double I = 0;
		for (iset.begin(); iset.dist < iset.size-1; ++iset){
			double unext = (use_log_densities)? exp(*(itu+1)) : *(itu+1);
			double unow  = (use_log_densities)? exp(*itu) : *itu;
			I += (*(itx+1)-*itx)*(w(*(itx+1), t)*unext + w(*itx, t)*unow);
		}
		return I*0.5;
	}
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


