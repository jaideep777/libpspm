#ifndef PSPM_PSPM_VECTOR_INSERT_H_
#define PSPM_PSPM_VECTOR_INSERT_H_

#include <vector>


// FIXME NOTE: ins_at must be sorted ascending. 
template <class T>
void vector_insert(std::vector<T> &x, std::vector<int>& ins_at, std::vector<T>& ins_vals){
	
	size_t n = x.size();
	x.resize(x.size()+ins_at.size());
	int at = ins_at.size()-1;
	int id = n-1;
	auto idx = x.begin()+x.size()-1;
	while(at >=0){
		//// below code is equivalent to:
		//while(id >= ins_at[at]){
		//    *idx = x[id];
		//    --idx; --id;
		//}
		if (ins_at[at]<=id){ // from first            from last+1     to last+1
			std::copy_backward(x.begin()+ins_at[at], x.begin()+id+1, idx+1);
			idx -= id-ins_at[at]+1;
			id  -= id-ins_at[at]+1;
		}
		*idx = ins_vals[at];
		--idx;
		--at;

		//for (auto xx :x) std::cout << xx << " ";
		//std::cout << std::endl;
	}

}

#endif
