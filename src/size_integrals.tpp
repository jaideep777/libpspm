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
	
	if (method == SOLVER_CM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the beginning, so x will be ascending
		double I = 0;
		double x_hi = *itx;
		double f_hi = w(*itx, t)*exp(*itu);
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		--itx; --itu;
		for (int i=0; i<spp.J-1; ++i){
			double x_lo = *itx;
			double f_lo = w(*itx,t)*exp(*itu);
			//if (xlow < 0.01) cout << "x/w/u/f = " << x_lo << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_lo << "\n";
			--itx; --itu;
			if (x_lo < xlow){
				double f = f_lo; //f_lo + (f_hi-f_lo)/(x_hi-x_lo)*(xlow - x_lo);  // FIXME: these should stop at the interpolating point
				double x = x_lo; //xlow;
				I += (x_hi-x) * (f_hi + f);
				break;
			}
			else{
				I += (x_hi - x_lo) * (f_hi + f_lo);
			}
			x_hi = x_lo;
			f_hi = f_lo;
		}
		
		//if (xsize() == 1 || f_hi > 0){
		//    double x_lo = xb;
		//    double g = mod->growthRate(xb, t);
		//    double u0 = (g>0)? u0_in*mod->establishmentProbability(t)/g  :  0; //FIXME: as of now, this does not work. To be fixed and discussed with Ake/Ulf
		//    double f_lo =  w(xb, t)*u0;
		//    I += (x_hi-x_lo)*(f_hi+f_lo);
		//}
		// for now, ignoring the case of single cohort - 0 will be returned, and 
		// that's probably okay.
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
		double * U = S.data();
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
		for (iset.begin(); !iset.end(); ++iset){
			I += (*(itx+1)-*itx)*(w(*(itx+1), t)*exp(*(itu+1)) + w(*itx, t)*exp(*itu));
		}
		return I*0.5;
	}
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


