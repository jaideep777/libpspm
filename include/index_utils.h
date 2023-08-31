#ifndef PSPM_INDEX_UTILS_H_
#define PSPM_INDEX_UTILS_H_

#include <vector>

// From tensorlib


	/// Convert 1D index to coordinates (Inverse of location())
	std::vector<int> index(int loc, std::vector<int> dim){
		int ndim = dim.size();
		std::vector<int> id(ndim);
		for (int i=ndim-1; i>=0; --i){
			int ix = loc % dim[i];
			loc = (loc-ix)/dim[i];
			id[i]=ix;
		}
		return id;
	}


    int location(std::vector <int> ix, std::vector<int> dim){
		int loc = 0;
		int ndim = dim.size();

        std::vector<int> offsets = std::vector<int>(ndim);
		int p = 1;
		for (int i=ndim-1; i>=0; --i){
			offsets[i] = p;
			p *= dim[i];
		}

		for (int i=ndim-1; i>=0; --i){
			loc += offsets[i]*ix[i];
		}
		return loc;
	}



#endif
