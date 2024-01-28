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

	if (method == SOLVER_EBT || method == SOLVER_IEBT){ 
		if (debug) std::cout << "~~~~ start step: " << current_time	<< " --> " << tstop << "\n";
		while(current_time < tstop){
			double tnext = std::fmin(tstop, t_next_cohort_insertion);

			if (debug) std::cout << std::setprecision(12) << "Stepping: " << current_time << " --> " << tnext << '\n';
			if (method == SOLVER_IEBT) stepTo_implicit(tnext, after_step, &Solver::stepU_iEBT);
			if (method == SOLVER_EBT)  stepTo_explicit(tnext, after_step, &Solver::calcRates_EBT);

			if (current_time >= t_next_cohort_insertion){ //} (fabs(current_time - t_next_cohort_insertion)<1e-12){
				if (debug) std::cout << std::setprecision(12) << "updating cohorts (" << current_time << ")\n";
				// update cohorts
				mergeCohorts_EBT();
				removeDeadCohorts_EBT();
				addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
				t_next_cohort_insertion += control.cohort_insertion_dt;
			}			

		}
	}

	if (method == SOLVER_CM || method == SOLVER_ICM){
		if (debug) std::cout << "~~~~ start step: " << current_time	<< " --> " << tstop << "\n";
		while(current_time < tstop){
			double tnext = std::fmin(tstop, t_next_cohort_insertion);

			if (debug) std::cout << std::setprecision(12) << "Stepping: " << current_time << " --> " << tnext << '\n';
			if (method == SOLVER_ICM) stepTo_implicit(tnext, after_step, &Solver::stepU_iCM);
			if (method == SOLVER_CM)  stepTo_explicit(tnext, after_step, &Solver::calcRates_CM);

			if (current_time >= t_next_cohort_insertion){ //} (fabs(current_time - t_next_cohort_insertion)<1e-12){
				if (debug) std::cout << std::setprecision(12) << "updating cohorts (" << current_time << ")\n";
				// update cohorts
				if (control.update_cohorts){
					addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
					removeCohort_CM();
				}
				t_next_cohort_insertion += control.cohort_insertion_dt;
			}			

		}
	}

	if (method == SOLVER_IFMU){
		stepTo_implicit(tstop, after_step, &Solver::stepU_iFMU);
	}

	if (method == SOLVER_FMU){	
		stepTo_explicit(tstop, after_step, &Solver::calcRates_FMU);
	}
	
	if (method == SOLVER_ABM){	
		stepTo_abm(tstop, after_step);
	}

	// std::cout << "Finished step to " <<std::endl;
}


template<typename AfterStepFunc>
void Solver::stepTo_explicit(double tstop, AfterStepFunc &afterStep, calcRatesPointer rates_func){
	auto derivs = [this, rates_func](double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt, void* params){
		if (debug) std::cout << "derivs()\n";
		copyStateToCohorts(S); 
		updateEnv(t, S, dSdt);
		(this->*rates_func)(t, S, dSdt);
	};
	
	// integrate
	odeStepper.step_to(tstop, current_time, state, derivs, afterStep); // rk4_stepsize is only used if method is "rk4"
}


// NOTE: This does a semi implicit update of system variables, as follows:
// This gives better reults that using only dSdt0 to step the s-state (tested with the Daphnia model) 
//   t0 -----------------> t+ -------------------> t1
//   X0,U0 --------------> X1,U1 ----------------> X1,U1
//   M0 -----------------> M1 -------------------> M1
//   S0 --------------------+-----------------+--> S1 = f((dSdt0 + dSdt+)/2)
//    +-- E0(X0,U0,S0) -->  +-- E+(X1,U1,S0) -^-->  +-- E1(X1,U1,S1)
//    +-- dSdt0(E0,S0) -->  +-- dSdt+(E+,S0) -^-->  +-- dSdt1(E0,S0)
// 
template<typename AfterStepFunc>
void Solver::stepTo_implicit(double tstop, AfterStepFunc &afterStep, stepUPointer step_u_func){
	while (current_time < tstop){
		double dt = std::fmin(control.ode_ifmu_stepsize, tstop-current_time);
		// std::cout << "   implicit step: " << current_time << " --> " << current_time + dt << "\n";
		
		//copyStateToCohorts(state.begin()); // not needed here because it is called by the odestepper below
		updateEnv(current_time, state.begin(), rates.begin()); // this computes E0, dSdt0
		std::vector<double> sys_rates_prev(rates.begin(), rates.begin()+n_statevars_system);  // save system rates (dSdt0)
		
		// use implicit stepper to advance u
		(this->*step_u_func)(current_time, state, rates, dt);
		// current_time += dt; // not needed here, as current time is advanced and state copied to cohorts by the ODE stepper below.
		// copyStateToCohorts(state.begin());   // copy updated X/U to cohorts 
		
		// Step accumulators before update of system vars because system vars will update env
		stepAccumulators(dt); // this will copyStateToCohorts before every ODE step (but not after the last step).
		copyStateToCohorts(state.begin()); // so copy once (because stepSystemVars needs to upate env)

		stepSystemVars(sys_rates_prev, dt);  // Computes E+, dSdt+, S1 (does not copy state at the end)
		afterStep(current_time, state.begin()); 
	}	
}


template<typename AfterStepFunc>
void Solver::stepTo_abm(double tstop, AfterStepFunc & afterStep){

	while (current_time < tstop){
		double dt = std::min(control.abm_stepsize, tstop-current_time);
		
		//copyStateToCohorts(state.begin()); // not needed here because it is called by the odestepper below
		updateEnv(current_time, state.begin(), rates.begin());
		std::vector<double> rates_prev(rates.begin(), rates.begin()+n_statevars_system);  // save system variable rates
		
		// use implicit stepper to advance u
		stepABM(current_time, dt);  // this will step all variables, including accumulators, and copies state to cohorts
		current_time += dt; 
		
		// step system vars
		stepSystemVars(rates_prev, dt);  // Computes E+, dSdt+, S1 (does not copy state at the end)

		// Need to explicitly call this because ODE solver is not used in ABM
		// Should cohorts be copied to state here? - done within stepABM() above
		afterStep(current_time, state.begin());
	}
}

