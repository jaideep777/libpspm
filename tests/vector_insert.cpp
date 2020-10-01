#include <iostream>


//template <class T>
//void vector_insert(std::vector<T> &x, std::vector<int>& ins_at, std::vector<T>& ins_vals){
	
//    size_t n = x.size();
//    x.resize(x.size()+ins_at.size());
//    int at = ins_at.size()-1;
//    int id = n-1;
//    auto idx = x.begin()+x.size()-1;
//    while(at >=0){
//        //// below code is equivalent to:
//        //while(id >= ins_at[at]){
//        //    *idx = x[id];
//        //    --idx; --id;
//        //}
//        if (ins_at[at]<id){ // from first            from last+1     to last+1
//            std::copy_backward(x.begin()+ins_at[at], x.begin()+id+1, idx+1);
//            idx -= id-ins_at[at]+1;
//            id  -= id-ins_at[at]+1;
//        }
//        *idx = ins_vals[at];
//        --idx;
//        --at;
	
//        //for (auto xx :x) std::cout << xx << " ";
//        //std::cout << std::endl;
	
//    }

//}

#include <vector_insert.h>

using namespace std;

int main(){

	
	vector <float> x = {1,2,3, 0.1, 0.2, 0.3,    6,7,8,9, .6,.7,.8,.9};
	
	vector <int> ins_at = {0,0,1,1, 3,3,4,4, 6, 6, 10, 10, 10, 14,14};
	vector <float> ins_vals = {.25,.5,1.25,1.5, .025,0.05,.125,.15, 4,5,10, .4,.5, 1.0,1.1}; 

	
	vector_insert(x, ins_at, ins_vals);


	float res[] = {.25, 0.5, 1, 1.25, 1.5, 2, 3, .025, 0.05, 0.1, .125, .15, 0.2, 0.3, 4, 5, 6, 7, 8, 9, 10, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1}; 
	for (int i=0; i<x.size(); ++i){
		if (x[i] != res[i]) return 1;
	}
	return 0;
	

}

