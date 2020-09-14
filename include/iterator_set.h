#ifndef PSPM_PSPM_ITERATOR_SET_H
#define PSPM_PSPM_ITERATOR_SET_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>

/**********************************************************
 *
 *   Iterator Set
 *
 *   This simple class provides a set of iterators to 
 *   the specified container.
 *
 *   It also provides for naming the iterators and 
 *   retrieving them by name.  
 *
 *   This is useful for traversing state arrays where 
 *   values of different variables are packed one after 
 *   the other in a custom layout.
 *
 * ********************************************************/


template <class Iterator>
class IteratorSet{
	public:
	Iterator first;
	std::vector<std::string> names;
	std::vector<Iterator> iters;
	std::vector<int> strides;
	int size;
	int dist = 0;
		
	public:
	// Create a set of equally spaced iterators starting at `first` and 
	// with names `varnames`, with spacing `size`
	IteratorSet(Iterator _first, std::vector<std::string> varnames, int _size, std::vector<int> spacing, std::vector<int> _strides){
		first = _first;
		size = _size;
		Iterator temp = first;
		names = varnames;
		strides = _strides;

		iters.push_back(first);
		for (int i=1; i<varnames.size(); ++i){
			iters.push_back(std::next(iters[i-1], spacing[i-1]));
		}
	}

	// Reset iterators
	void begin(){
		for (int i=0; i<names.size(); ++i){
			iters[i] -= dist*strides[i];
		}
		dist = 0;
	}

	// Check whether the iterator set has traversed till the end
	bool end(){
		return dist >= size; //iters[0] == next(first,strides[0]*size);
	}

	// add a new iterator to the set
	void push_back(std::string name, Iterator it, int stride){
		begin();
		names.push_back(name);
		iters.push_back(it);
		strides.push_back(stride);
	}
	
	// get all iterators in a vector. Returns a const reference to 
	// the internal vector of iterators so that it remains valid 
	// even after ++ is called on the original
	const std::vector<Iterator>& get(){
		return iters;
	}

	// increment all iterators
	IteratorSet& operator++(){
		for (int i=0; i<iters.size(); ++i) iters[i] += strides[i];
		++dist;
		return *this;
	}

	// Get a specific iterator by name
    const Iterator& get(std::string name){
		int id;
		for (id=0; id<names.size(); ++id) if (names[id] == name) break;
        return iters[id];    
    }

	// Get a specific iterator index by name
    size_t getIndex(std::string name){
		int id;
		for (id=0; id<names.size(); ++id) if (names[id] == name) break;
        return id;    
    }
	
	// Print formatted iterator names and values
	void printHeader(int w = 11){
		for (auto n : names) std::cout << std::setw(w) << n << " "; 
		std::cout << "\n";
	}
	
	void printLine(int w = 11){
		for (auto& it : iters) std::cout << std::setw(w) << *it << " ";
		std::cout << "\n";
	}

	void print(unsigned int n = -1){
		printHeader();
		for (begin(); !end() && n>0; ++(*this), --n){
			printLine();
		}
	}

	void printInfo(){
		std::cout << "--- Iterator Set ---\n";
		for (int i=0; i<names.size(); ++i){
			std::cout << names[i] << "\t" << strides[i] << " :\t" << *iters[i] << "\n";
		}
		std::cout << "--------------------\n";

	}
};


#endif



