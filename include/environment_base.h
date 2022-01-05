#ifndef  PSPM_PSPM_ENVIRONMENT_H_
#define  PSPM_PSPM_ENVIRONMENT_H_

class Solver;

class EnvironmentBase{
	public:
	virtual void computeEnv(double t, Solver * sol) = 0;

};


#endif