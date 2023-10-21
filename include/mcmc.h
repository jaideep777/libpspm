#ifndef PSPM_MCMC_H_
#define PSPM_MCMC_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>

class Chain{
	public:
	int dim;
	std::vector<std::vector<double>> samples; // each sample is an nD vector
	std::mt19937 rng;

	public:
	Chain(const std::vector<double>& start) : rng(std::random_device()()){
		dim = start.size();
		samples.push_back(start);
	}

	std::vector<double> proposal(const std::vector<double>& current, const std::vector<double>& sd) {
		std::vector<double> proposed(dim);
		for (int i=0; i<dim; ++i) {
			std::normal_distribution<double> proposal_dist(current[i], sd[i]);
			proposed[i] = proposal_dist(rng);
		}
		return proposed;
	}

	// Func f must be a function that takes a vector x as input and returns a double
	template <class Func>
	std::vector<double> sample(Func target, const std::vector<double>& sd){
		// Generate a proposal from the multivariate proposal distribution
		std::vector<double> x_current = *samples.rbegin();
		std::vector<double> x_proposed = proposal(x_current, sd);

		double rand = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
		if (rand < target(x_proposed) / target(x_current)){
			samples.push_back(x_proposed);
		}
		else{
			samples.push_back(x_current);
		}
		return *samples.rbegin();	
	}

	void printSamples(std::ostream &fout){
		for (int i=0; i<samples.size(); ++i){
			fout << i << '\t';
			for (auto c : samples[i]) fout << c << '\t';
			fout << '\n';
		}
	}
};

// This class should initialize and manage chains, 
// deal with convergence, and provide samples.
class MCMCSampler {

};



#endif

