#include <iostream>
#include "../include/mcmc.h"
using namespace std;

// compile: g++ -o 1 mcmc_test.cpp

// Test with R:
// ------------------
// library(tidyverse)
// dat = read.delim("/home/jjoshi/codes/libpspm/tests/mcmc_test.txt", header=F)
// dat %>% ggplot(aes(x=V2, y=V3, col=V1)) + 
//   geom_point(alpha=0.5) + 
//   # geom_line(alpha=0.1) +
//   geom_point(x=1,y=2,col="red")+
//   geom_point(x=5,y=2,col="red")+
//   geom_point(x=1,y=10,col="red")+
//   geom_point(x=5,y=10,col="red")

// hist(dat$V2[-(1:200)], breaks=50)
// hist(dat$V3[-(1:200)], breaks=50)


int main() {

	auto targetDistribution = [](const std::vector<double>& x) {
        std::vector<double> means = {1,2};
		double result = 1.0;
        for (int i=0; i<x.size(); ++i) {
            result *= std::exp(-(x[i]-means[i])* (x[i]-means[i])) + std::exp(-(x[i]-5*means[i])* (x[i]-5*means[i])); // Example: Multivariate Normal Distribution
        }
        return result;
    };

	Chain C({15,0});
	for (int s=0; s<10000; ++s){
		C.sample(targetDistribution, {.5,.5});
	}

	Chain C2({15,15});
	for (int s=0; s<10000; ++s){
		C2.sample(targetDistribution, {.5,.5});
	}

	std::ofstream fout("mcmc_test.txt");
	C.printSamples(fout);
	C2.printSamples(fout);
	fout.close();

	return 0;
}
