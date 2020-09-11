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
 *   This simple class provides a set of equally spaced
 *   iterators to the specified container.
 *
 *   It also provides for naming the iterators and 
 *   retrieving them by name.  
 *
 *   This is useful for traversing state arrays where 
 *   values of different variables are packed one after 
 *   the other. 
 *
 * ********************************************************/


template <class Iterator>
class IteratorSet{
	private:
	Iterator first;
	std::vector<std::string> names;
	std::vector<Iterator> iters;
	std::map<std::string, size_t> names_map;
	int size;
		
	public:
	// Create a set of equally spaced iterators starting at `first` and 
	// with names `varnames`, with spacing `size`
	IteratorSet(Iterator _first, std::vector<std::string> varnames, int stride){
		first = _first;
		size = stride;
		for (int i=0; i<varnames.size(); ++i){
			names.push_back(varnames[i]);
			iters.push_back(std::next(first, stride*i));
			names_map[varnames[i]] = i;
		}
	}

	// Reset iterators
	void begin(){
		for (int i=0; i<names.size(); ++i){
			iters[i] = std::next(first, size*i);
		}
	}

	// Check whether the iterator set has traversed till the end
	bool end(){
		return iters[0] == next(first,size);
	}

	// add a new iterator to the set
	void push(std::string name, Iterator it){
		names.push_back(name);
		iters.push_back(it);
		names_map[name] = iters.size()-1;
	}
	
	// get all iterators in a vector. Returns a const reference to 
	// the internal vector of iterators so that it remains valid 
	// even after ++ is called on the original
	const std::vector<Iterator>& get(){
		return iters;
	}

	// increment all iterators
	IteratorSet& operator++(){
		for (auto& it : iters) ++it;
		return *this;
	}

	// Get a specific iterator by name
    const Iterator& get(std::string name){
        return iters[names_map[name]];    
    }

	// Print formatted iterator names and values
	void printHeader(int w = 11){
		for (auto n : names) std::cout << std::setw(w) << n << " "; 
		std::cout << "\n";
	}
	
	void print(int w = 11){
		for (auto& it : iters) std::cout << std::setw(w) << *it << " ";
		std::cout << "\n";
	}



};


#endif



