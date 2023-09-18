#ifndef PSPM_INDEX_UTILS_H_
#define PSPM_INDEX_UTILS_H_

#include <vector>

// From tensorlib

namespace id_utils{


	/// Convert 1D index to coordinates (Inverse of location())
	std::vector<int> index(int loc, const std::vector<int>& dim){
		int ndim = dim.size();
		std::vector<int> id(ndim);
		for (int i=ndim-1; i>=0; --i){
			int ix = loc % dim[i];
			loc = (loc-ix)/dim[i];
			id[i]=ix;
		}
		return id;
	}


    int location(std::vector <int> ix, const std::vector<int>& dim){
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

//       |
//  44   |
//  33   |------x  : id = [1,2]    
//  22   |      .    coords = [[10,20*,30,...]       coords[1]
//  11   |      .              [11,22,33*,44,...]]   coords[2]
//       |--------------------------------
//         10   20   30   40
std::vector<double> coord_value(const std::vector <int>& id, const std::vector<std::vector<double>>& coords){
	std::vector<double> vals(id.size());
	for (int i=0; i<id.size(); ++i){
		vals[i] = coords[i][id[i]];
	}
	return vals;
}

}


#endif
