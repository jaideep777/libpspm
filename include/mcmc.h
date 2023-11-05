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

	public:
	std::vector<Chain> chainList;
	int num_Chains;
	int dims; 
	int burn_in; 
	int thinning; 
	const std::vector<double> x_min;
	const std::vector<double> x_max; 
	const std::vector<double> sd; 
	std::vector<std::vector<double>> merged_chains;

	MCMCSampler(const std::vector<double> _x_min, 
				const std::vector<double> _x_max, 
				const std::vector<double>& _sd,
				int _nChains, 
				int _burn_in, 
				int _thinning): 
				x_min(_x_min), 
				x_max(_x_max),
				sd(_sd){

		std::default_random_engine generator;

		num_Chains = _nChains;
		burn_in = _burn_in;
		thinning = _thinning;
		dims = x_min.size();

		for(int i = 0; i < num_Chains; ++i){
			std::vector<double> start;
			for(int k = 0; k < x_min.size(); k++){
				std::uniform_real_distribution<> random_start(x_min[k],x_max[k]);
				start.push_back(random_start(generator));
			}
			chainList.push_back(Chain(start));
			merged_chains.push_back(start);
		}
	}

	template <class Func>
	void run_chains(Func target, int chainLength){
		for(int i=0; i < chainLength; ++i){
			for (auto & chain : chainList) 
  			{
				while (!accept(chain.sample(target,sd)))
				{
					/* do nothing */
				}
				merged_chains.push_back(chain.samples.back());
  			}
		}
	}

	std::vector<double> gelman_rubin_test(){
		
		std::vector<std::vector<double>> posterior_mean_chain(dims, std::vector<double> (num_Chains, 0.0));
		std::vector<double> posterior_mean(dims);
		std::vector<std::vector<double>> intra_chain_variance_chain(dims, std::vector<double> (num_Chains, 0.0));

		int num_elements = chainList[0].samples.size() - burn_in;

		std::vector<double> B(dims, 0); // B = how the individual means vary around the joint mean
		std::vector<double> W(dims, 0); // W = averaged variances of the chains

		std::vector<double> V(dims, 0); // V = true variance
		std::vector<double> R(dims, 0); // R = test for convergence

		for (int k = 0; k < dims; ++k)
		{
			
			for (int j = 0; j < num_Chains; ++j)
			{
				for (int i = burn_in; i < chainList[0].samples.size(); ++i)
				{
					posterior_mean_chain[k][j] += chainList[j].samples[i][k] / num_elements;
				}

				for (int i = burn_in; i < chainList[0].samples.size(); ++i)
				{
					intra_chain_variance_chain[k][j] += pow((chainList[j].samples[i][k] - posterior_mean_chain[k][j]),2) / (num_elements-1);
				}

				posterior_mean[k] += posterior_mean_chain[k][j] / num_Chains;
				W[k] += intra_chain_variance_chain[k][j] / num_Chains;
			}

			for (int j = 0; j < num_Chains; ++j)
			{
				B[k] += pow((posterior_mean_chain[k][j] - posterior_mean[k]), 2) * num_elements / (num_Chains-1);
			}

			V[k] = 1.0 * (num_elements - 1)/num_elements * W[k] + (num_Chains + 1)/(num_Chains * num_elements) * B[k];
			R[k] = sqrt(V[k]/W[k]); // should ~~ 1 - need a tolerance measure
		}

		// test Rk - for now just returning R because... well 
		return R;
	}

	std::vector<std::vector<double>> sample(int nSamples){
		std::vector<std::vector<double>>::const_iterator first = merged_chains.end() - nSamples + 1;
		std::vector<std::vector<double>>::const_iterator last = merged_chains.end();
		std::vector<std::vector<double>> newVec(first, last);
		return newVec;
	}

	bool accept(std::vector<double> element){
		bool out = true;
		for(int i = 0; i < dims ; ++i){
			if(element[i] < x_min[i] || element[i] > x_max[i]){
				return false;
			}
		}
		return out;
	}

};



#endif

