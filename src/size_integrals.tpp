template<class Model>
template<typename wFunc>
double Solver<Model>::integrate_wudx_above(wFunc w, double t, double xlow, vector<double>&S){
	//cout << " | " <<  t << " " << mod->evalEnv(0,t) << " ";
	//if (method == SOLVER_FMU){
	//}
	//else if (method == SOLVER_EBT){
	//}
	//cout << "Begin integrate: xsize = " << xsize() << "(" << S[0] << ", " << S[xsize()-1] << "), xlow = " << xlow << endl;
	
	if (method == SOLVER_CM){
		// integrate using trapezoidal rule 
		// Note, new cohorts are inserted at the beginning, so x will be ascending
		auto itx = S.begin() + xsize()-1;
		auto itu = S.begin() + 2*xsize()-1;
		double I = 0;
		double x_hi = *itx;
		double f_hi = w(*itx, t)*exp(*itu);
		//if (xlow < 0.01) cout << "x/w/u/f = " << x_hi << " " <<  w(*itx,t) <<  " " << exp(*itu)  << " " << f_hi << "\n";
		--itx; --itu;
		for (int i=0; i<xsize()-1; ++i){
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






template<class Model>
template<typename wFunc>
double Solver<Model>::integrate_x(wFunc w, double t, vector<double>&S, int power){
	//cout << " | " <<  t << " " << mod->evalEnv(0,t) << " ";
	if (method == SOLVER_FMU){
		// integrate using midpoint quadrature rule
		double I=0;
		double * U = S.data();
		for (unsigned int i=0; i<X.size(); ++i){
			I += h[i]*w(X[i], t)*pow(U[i], power);  // TODO: Replace with std::transform after profiling
		}
		return I;
	}
	else if (method == SOLVER_EBT){
		// integrate using EBT rule (sum over cohorts)
		double   pi0  =  S[0];
		double * xint = &S[1];
		double   N0   =  S[xsize()];
		double * Nint = &S[xsize()+1];
		
		double x0 = xb + pi0/(N0+1e-12); 
		
		double I = w(x0, t)*N0;
		for (int i=0; i<xsize()-1; ++i) I += w(xint[i], t)*Nint[i];
		
		return I;
	}
	else if (method == SOLVER_CM){
		// integrate using trapezoidal rule TODO: Modify to avoid double computation of w(x)
		double * px = &S[0];
		double * pu = &S[xsize()];
		double I = 0;
		for (int i=0; i<xsize()-1; ++i){
			I += (px[i+1]-px[i])*(w(px[i+1], t)*exp(pu[i+1]) + w(px[i], t)*exp(pu[i]));
		}
		return I*0.5;
	}
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


