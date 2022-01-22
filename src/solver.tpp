// current_time is updated by the ODE solver at every (internal) step
template<typename AfterStepFunc>
void Solver::step_to(double tstop, AfterStepFunc &afterStep_user){
	// do nothing if tstop is <= current_time
	if (tstop <= current_time) return;
	
	auto after_step = [this, afterStep_user](double t, std::vector<double>::iterator S){
		//cout << "After step: t = " << t << "\n";
		copyStateToCohorts(S);
		//for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	// FIXME: Check if this is needed
		//for (auto spp : species_vec) spp->afterStep(t, env);
		afterStep_user(t);
	};

	
	if (method == SOLVER_FMU){	
		auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
			copyStateToCohorts(S);
			env->computeEnv(t, this);
			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			this->calcRates_FMU(t, S, dSdt);
		};
		
		odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
//		copyStateToCohorts(state.begin());
	}
	
	
	if (method == SOLVER_IFMU){	
		while (current_time < tstop){
			double dt = std::min(control.ode_ifmu_stepsize, tstop-current_time);
			
			//copyStateToCohorts(state.begin());
			env->computeEnv(current_time, this);
			
			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,current_time);	
			
			// use implicit stepper to advance u
			stepU_iFMU(current_time, state, rates, dt);
			// current_time += dt; // not needed here, as current time is advanced by the ODE stepper below.

			// use the ODE-stepper for other state variables
			auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
				copyStateToCohorts(S);
				// precompute and env computation is not needed here, because it depends on x and u, which are not updated by the solver.
				calcRates_iFMU(t,S,dSdt);
			};
			// this step below will do afterstep. FIXME: But what if there are no extra state variables? 
			odeStepper.step_to(current_time+dt, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
	
//			copyStateToCohorts(state.begin());
		}

	}
	
	//if (method == SOLVER_MMU){
	//}
	
	
	if (method == SOLVER_EBT){
		auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
			// copy state std::vector to cohorts
			copyStateToCohorts(S);
			
			// compute environment
			env->computeEnv(t, this);

			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			// get rates
			calcRates_EBT(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
		
		// after the last ODE step, the state vector is updated but cohorts still hold an intenal ODE state (y+k5*h etc).
		// normally, this will be no problem since state will be copied to cohorts in the next rates call. 
		// But since add/remove cohort below will rewrite the state from cohorts, the updated state vector will be lost
		// rewrite the cohorts now to avoid this.
//		copyStateToCohorts(state.begin());
		
		// update cohorts
		removeDeadCohorts_EBT();
		addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
	}
	
	
	if (method == SOLVER_CM){
		auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
			// copy state vector to cohorts
			copyStateToCohorts(S);

			// update u0 (u of boundary cohort)
			for (auto s : species_vec) s->get_u0(t, env);
			
			// compute environment
			env->computeEnv(t, this);

			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			// get rates
			calcRates_CM(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
		
		// after the last ODE step, the state vector is updated but cohorts still hold an intenal ODE state (y+k5*h etc).
		// normally, this will be no problem since state will be copied to cohorts in the next rates call. 
		// But since add/remove cohort below will rewrite the state from cohorts, the updated state vector will be lost
		// rewrite the cohorts now to avoid this.
		// Note: This is likely no longer required because afterStep() does a copy
		copyStateToCohorts(state.begin());

		// update cohorts
		if (control.update_cohorts){
			addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
			removeCohort_CM();
		}
		//env->computeEnv(current_time, this); // is required here IF rescaleEnv is used in derivs
	}
}

