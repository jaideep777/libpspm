#include <cassert>

// ~~~~~~~~~~~ FMU Solver ~~~~~~~~~~~
double phi(double r){
	return max(max(0.0,min(2*r,1.0)),min(r,2.0));
}


template<class Model, class Environment>
void Solver<Model,Environment>::calcRates_FMU(double t, vector<double> &S, vector<double> &dSdt){

	for (int species_id=0; species_id < species_vec.size(); ++species_id){
		Species<Model> &spp = species_vec[species_id];

		//cout << spp.start_index << endl;
		double *U = &S[spp.start_index];		// offset for species
		double *dUdt = &dSdt[spp.start_index];	// offset for species
		int J = spp.J;			// xsize of species.
		vector<double> &x = spp.x;
		vector<double> &X = spp.X;
		vector<double> &h = spp.h;

		// for something like Plant, we will need to calc mortality rate and extra rates immediately after growthrate. 
		// need extra arrays to store, or ensure order is g[i], g[i+1/2], g[i+1]
		vector <double> growthArray(J+1);
		for (int i=0; i<J+1; ++i) growthArray[i] = spp.mod->growthRate(x[i], t, env);

	//	#define growth(i) growthRate(x[i], mod->evalEnv(x[i],t))
		#define growth(i) growthArray[i]

		vector <double> u(J+1);
		
		// i=0
		if (spp.birth_flux_in < 0){	
			double birthFlux = integrate_x([this, species_id](double z, double t){return species_vec[species_id].mod->birthRate(z,t, env);}, t, S, species_id);
			u[0] = birthFlux/(growth(0)+1e-12); // Q: is this correct? or g(X0,env)? - A: It is g(xb,env) - growth() indexes x, so (0) is xb. Hence correct. 
		}
		else{
			u[0] = spp.birth_flux_in;  // FIXME: for now, using birth_flux_in as u0_in for test_model
		}

		// i=1 (calc u1 assuming linear u(x) in first interval)
		u[1] = 2*U[0]-u[0];  // NOTE: for g(x) < 0 this can be calculated with upwind scheme 
		
		for (int i=2; i<J-1; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			if(growth(i) >=0){
				double rMinus = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i-1]-U[i-2]+1e-12)/(x[i-1]-x[i-2]));
				u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
			}   
			else{
				double rPlus  = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i+1]-U[i]+1e-12)/(x[i+1]-x[i]));
				u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]); 
			}
		}
		
		u[J-1] = 2*U[J-2] - u[J-2];	// NOTE: for g(x) > 0 This can be calc with upwind scheme
		u[J] = 2*U[J-1] - u[J-1];

		auto is = spp.get_iterators(S);
		auto ir = spp.get_iterators(dSdt);
		auto& itx = is.get("X");
		auto& itu = is.get("u");
		auto& itse = is.get((spp.varnames_extra.size()>0)? spp.varnames_extra[0] : "X");	// dummy init //FIXME: commented line below does not work. Explore why
		//auto& itse = (varnames_extra.size()>0)? is.get(varnames_extra[0]) : S.begin(); 
		
		auto& itdu = ir.get("u");
		auto& itre = ir.get((spp.varnames_extra.size()>0)? spp.varnames_extra[0] : "X");	// dummy init
		//auto& itre = (varnames_extra.size()>0)? ir.get(varnames_extra[0]) : dSdt.begin(); 
		
		is.begin(), ir.begin();
		for (int i=0; !is.end(); ++i, ++is, ++ir){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			dUdt[i] = -spp.mod->mortalityRate(X[i], t, env)*U[i] - (growth(i+1)*u[i+1] - growth(i)*u[i])/h[i];
			
			if (spp.varnames_extra.size() > 0){
				auto it_returned = spp.mod->calcRates_extra(t, *itx, itse, itre);
				assert(distance(itre, it_returned) == spp.varnames_extra.size());
			}
			
		}

	}
}

	

