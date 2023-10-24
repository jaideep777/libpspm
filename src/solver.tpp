/// @param tstop   The time until which the ODE solver should be stepped. 
///                The stepper will stop exactly at tstop, for which the final step size is truncated if necessary.
/// @param afterStep_user  A function of the form `f(double t)` to be called after every _successful_ ODE step. 
///
/// \image html ode_flow.png width=700cm 
/// @note 1. `current_time` is updated by the ODE solver at every (internal) step
///
/// @note 2. After the last ODE step, the state vector is updated but cohorts still hold an intenal ODE state (y+k5*h etc).
/// normally, this will not be a problem because state will be copied to cohorts before rates call of the next timestep. 
/// But since add/remove cohort after step_to will rewrite the state from cohorts, the updated state vector will be lost.
/// To avoid this, we should ensure that the state is copied to cohorts after every successful ODE step (or at least
/// after completion of `step_to`).
/// This is achieved by the afterStep function
template<typename AfterStepFunc>
void Solver::step_to(double tstop, AfterStepFunc &afterStep_user){
	// do nothing if tstop is <= current_time
	// std::cout << "step to: current time: " << current_time << "\tt_stop: " << tstop << std::endl;

	if (tstop <= current_time) return;
	
	// std::cout << "Running step to function " << std::endl;

	auto after_step = [this, afterStep_user](double t, std::vector<double>::iterator S){
		if (debug) std::cout << "After step: t = " << t << "\n";
		copyStateToCohorts(S);
		afterStep_user(t);
	};

	
	if (method == SOLVER_IEBT){//} || method == SOLVER_IFMU || method == SOLVER_ICM){	
		while (current_time < tstop){
			double dt = std::min(control.ode_ifmu_stepsize, tstop-current_time);
			
			//copyStateToCohorts(state.begin()); // not needed here because it is called by the odestepper below
			updateEnv(current_time, state.begin(), rates.begin());
			std::vector<double> rates_prev(rates.begin(), rates.begin()+n_statevars_system);  // save system variable rates
			
			// use implicit stepper to advance u
			if      (method == SOLVER_IEBT) stepU_iEBT(current_time, state, rates, dt);
			// else if (method == SOLVER_IFMU) stepU_iFMU(current_time, state, rates, dt);
			// else if (method == SOLVER_ICM)  stepU_iCM(current_time, state, rates, dt);
			else     throw std::runtime_error("step_to(): Invalid solver method"); // TODO: This will actually never be reached so we can refactor to remove 
			// current_time += dt; // not needed here, as current time is advanced by the ODE stepper below.
			copyStateToCohorts(state.begin());   // copy updated X/U to cohorts 
			
			if (n_statevars_system > 0){ // FIXME: maybe keep this simple and just use the env and rates at be beginning of the timestep?
				// copyStateToCohorts(state.begin());   // not needed here as it's just done above
				updateEnv(current_time, state.begin(), rates.begin());  // FIXME: recompute env with updated u <-- is this necessary or even correct??
				// .FIXME: use fully implicit stepper here?
				for (int i=0; i<n_statevars_system; ++i){
					state[i] += (rates_prev[i]+rates[i])/2*dt;  // use average of old and updated rates for stepping system vars
				}
			}

			// FIXME: Move this to before system vars update, otherwise, this will see updated Env if sys vars are defined
			// FIXME: use an if condition to do this only if extra i-state is requested
			// use the ODE-stepper for other state variables
			// this stepper is called even if there are no extra state variables, so copyStateToCohorts is accomplished here
			// FIXME: Keep this simple and just take 1 Euler step?
			auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
				copyStateToCohorts(S);
				// precompute and env computation is not needed here, because it depends on x and u, which are not updated by the solver.
				calcOdeRatesImplicit(t,S,dSdt);
			};
			// this step below will do afterstep. -- Now doing explicitly, because when timestep is specified, we want afterstep only after the timestep is complete
			// Q: But what if there are no extra state variables? Still works, because ODE stepper will still step X/U with zero rates
			auto afterStep_dummy = [](double t, std::vector<double>::iterator S){};
			odeStepper.step_to(current_time+dt, current_time, state, derivs, afterStep_dummy); // rk4_stepsize is only used if method is "rk4"
	
			after_step(current_time, state.begin());
		}

		if (method == SOLVER_IEBT){
			// update cohorts
			mergeCohorts_EBT();
			removeDeadCohorts_EBT();
			addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
		}
		// if (method == SOLVER_ICM){
		// 	// update cohorts
		// 	if (control.update_cohorts){
		// 		addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
		// 		removeCohort_CM();
		// 	}
		// }


	}

	// if (method == SOLVER_FMU){	
	// 	auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
	// 		copyStateToCohorts(S); 
	// 		updateEnv(t, S, dSdt);
	// 		calcRates_FMU(t, S, dSdt);
	// 	};
		
	// 	// integrate
	// 	odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
	// }
	
	if (method == SOLVER_EBT){
		auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
			copyStateToCohorts(S);
			updateEnv(t, S, dSdt);
			calcRates_EBT(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
		
		// update cohorts
		mergeCohorts_EBT();
		removeDeadCohorts_EBT();
		addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
	}
	
	// if (method == SOLVER_EBTN){
	// 	auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
	// 		// std::cout << "in derivs: " << "\tstarting" << std::endl;
	// 		copyStateToCohorts(S);
	// 		// std::cout << "in derivs: " << "\tfinished copying state to cohorts" << std::endl;
	// 		updateEnv(t, S, dSdt);
	// 		// std::cout << "in derivs: " << "\tfinished updating environment" << std::endl;
	// 		calcRates_EBTN(t, S, dSdt);
	// 		// std::cout << "in derivs: " << "\tfinished calculating rates EBTN" << std::endl;
	// 	};
		
	// 	// integrate 
	// 	odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
		
	// 	// update cohorts
	// 	// mergeCohorts_EBTN();
	// 	removeDeadCohorts_EBTN();
	// 	addCohort_EBTN();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
	// }
	
	
	// if (method == SOLVER_CM){
	// 	auto derivs = [this](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
	// 		if (debug) std::cout << "derivs()\n";
	// 		copyStateToCohorts(S);  // this triggers precompute
	// 		updateEnv(t, S, dSdt);
	// 		calcRates_CM(t, S, dSdt);
	// 	};
		
	// 	// integrate 
	// 	odeStepper.step_to(tstop, current_time, state, derivs, after_step); // rk4_stepsize is only used if method is "rk4"
		
	// 	// update cohorts
	// 	if (control.update_cohorts){
	// 		addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
	// 		removeCohort_CM();
	// 	}
	// 	//env->computeEnv(current_time, this); // is required here IF rescaleEnv is used in derivs
	// }
	
	// if (method == SOLVER_ABM){	
	// 	while (current_time < tstop){
	// 		double dt = std::min(control.abm_stepsize, tstop-current_time);
			
	// 		//copyStateToCohorts(state.begin()); // not needed here because it is called by the odestepper below
	// 		updateEnv(current_time, state.begin(), rates.begin());
	// 		std::vector<double> rates_prev(rates.begin(), rates.begin()+n_statevars_system);  // save system variable rates
			
	// 		// use implicit stepper to advance u
	// 		stepABM(current_time, dt);  // this will step all variables, including extra_istate
	// 		current_time += dt; 
			
	// 		// step system vars
	// 		if (n_statevars_system > 0){
	// 			updateEnv(current_time, state.begin(), rates.begin());  // recompute env with updated u
	// 			// .FIXME: use fully implicit stepper here?
	// 			for (int i=0; i<n_statevars_system; ++i){
	// 				state[i] += (rates_prev[i]+rates[i])/2*dt;  // use average of old and updated rates for stepping system vars
	// 			}
	// 		}

	// 		// Need to explicitly call this because ODE solver is not used in ABM
	// 		// Should cohorts be copied to state here? - done within stepABM() above
	// 		after_step(current_time, state.begin());
	// 	}

	// }

	// std::cout << "Finished step to " <<std::endl;

}

