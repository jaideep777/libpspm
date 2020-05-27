


template <class _Model>
void Solver::initialize(_Model mod){
	// state vector was initialized to 0 in Constrctor. Set non-zero elements here
	switch (method){
		case SOLVER_FMU:	
			for (size_t i=0; i<J; ++i)  state[i] = mod.calcIC(X[i]);
			break;

		case SOLVER_MMU:
			for (size_t i=0; i<J; ++i)  state[J+1 + i] = mod.calcIC(X[i]);
			//for (size_t i=0; i<J+1; ++i) uprev[i] = calcIC(x[i]);
			break;

		case SOLVER_CM:
			for (size_t i=0; i<J+1; ++i)  state[J+1 + i] = mod.calcIC(x[i]);
			break;

		case SOLVER_EBT:
			for (size_t i=0; i<J; ++i)  state[J+1 + 1+i] = mod.calcIC(X[1+i]);	// state[J+1+0]=0 (N0)
			break;
	}
}



