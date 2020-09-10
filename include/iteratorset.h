#ifndef PSPM_PSPM_ITERATOR_SET_H
#define PSPM_PSPM_ITERATOR_SET_H

#include <vector>
#include <string>
#include <map>
#include <iostream>

template <class Iterator>
class IteratorSet{
	public:
	//std::vector<Iterator> iterators;
	std::map<std::string, Iterator> iters;
	
	public:
	IteratorSet(Iterator first, int num_vars, std::vector<std::string> names, int stride){
		for (int i=0; i<num_vars; ++i){
			iters[names[i]] = std::next(first, stride*i);
		}
	}
	
	std::map<std::string, Iterator> getIterators_map(){
		return iters;
	}
	
	IteratorSet& operator++(){
		for (auto& pairr : iters) ++pairr.second;
		return *this;
	}

    const Iterator& get(std::string name){
        return iters[name];    
    }
    
};


#endif



//int main(){
//    std::vector<double> state = {0,0.1,0.2,0.3, 1,1.1,1.2,1.3, 2,2.1,2.2,2.3, 3,3.1,3.2,3.3};
//    std::vector<std::string> names = {"0", "1", "2", "3"};
//    IteratorSet<std::vector<double>::iterator> iset(state.begin(), 4, names, 4);
//    
//    auto& it_0 = iset.get("0");
//    auto& it_1 = iset.get("1");
//    auto& it_2 = iset.get("2");
//    auto& it_3 = iset.get("3");
//    for (int i=0; i<4; ++i){
//        std::cout << *it_0 << " " << *it_1 << " " << *it_2 << " " << *it_3 << "\n";
//        ++iset;
//    }
//    
//}

