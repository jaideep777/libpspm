#ifndef MATH_TENSOR_H_
#define MATH_TENSOR_H_

#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>

#include <numeric>


/**
 Tensor. Multidimensional Array
 
 A tensor with dimensions
 ```
 { N2, N1, N0 } 
 ```
 is as follows.
 \image html tensorlib_tensor.png width=400cm 
 
 Elements of this tensor are accessed with a vector of indices as
 ```
 { i2, i1, i0 }
 ```
 i0 is the lowest dimension, i.e., (elements along i0 are stored consequtively in memory.
  This order of indices is chosen rather than {i0, i1, i2} to allow matrices 
  (two dimensional tensors) to be refered as {irow, icolumn}. 

 Thus the indices of each element of a 2x3x5 Tensor are 
 ``` 
   loc: index 
   ---:------
     0: 0 0 0 
     1: 0 0 1 
     2: 0 0 2 
     3: 0 0 3 
     4: 0 0 4 
     5: 0 1 0 
     6: 0 1 1 
          :   
    25: 1 2 0 
    26: 1 2 1 
    27: 1 2 2 
    28: 1 2 3 
    29: 1 2 4  
 ```   
 
 */


template <class T>
class Tensor{
	private:
	std::vector<int> offsets;

	protected:
	int nelem;
	
	public:
	std::vector<int> dim;
	std::vector<T> vec;

	/// Create a tensor with specified dimensions.
	/// This function also allocates space for the tensor, and calculates the offsets used for indexing.
	Tensor(std::vector<int> _dim){
		dim = _dim;
		nelem = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<int>());
		vec.resize(nelem);

		int ndim = dim.size();
		offsets.resize(ndim,0);
		int p = 1;
		for (int i=ndim-1; i>=0; --i){
			offsets[i] = p;
			p *= dim[i];
		}
		
	}


	/// Print the tensor.
	/// If vals is true, then values are also printed. Otherwise, only metadata is printed.
	void print(bool vals = true){
	    std::cout << "Tensor:\n";
	    std::cout << "   dims = "; for (auto d : dim) std::cout << d << " "; std::cout << "\n";
	    std::cout << "   offs = "; for (auto d : offsets) std::cout << d << " "; std::cout << "\n";
		if (vals){
			std::cout << "   vals = \n      "; std::cout.flush();
			for (int i=0; i<nelem; ++i){
				std::cout << vec[i] << " "; 
				bool flag = true;
				for (int axis=dim.size()-1; axis>0; --axis){
					flag = flag && (index(i)[axis] == dim[axis]-1);
					if (flag) std::cout << "\n      ";
				}
			}
		}
		std::cout << "\n";
	}
	
//	TODO: 
//	This function can be private
	/// Convert coordinates (specified as a vector of indices) to 1D index where the value resides in the underlying vector.
	int location(std::vector <int> ix){
		int loc = 0;
		int ndim = dim.size();
		for (int i=ndim-1; i>=0; --i){
			loc += offsets[i]*ix[i];
		}
		return loc;
	}

	template<class... ARGS>
	int location(ARGS... ids){
		return location({ids...});
	}


	/// @brief Get the value at coordinates specified as a comma separated list. 
	///        The order of coordinates is \f$\{i_n, i_{n-1}, ..., i_0\}\f$
	template<class... ARGS>
	double& operator() (ARGS... ids){
		return vec[location({ids...})];
	}

	/// @brief Get the value at coordinates specified as an integer vector. 
	///        The order of coordinates is \f$\{i_n, i_{n-1}, ..., i_0\}\f$.
	double& operator() (std::vector<int> ix){
		return vec[location(ix)];
	}


	/// Convert 1D index to coordinates (Inverse of location())
	std::vector<int> index(int loc){
		int ndim = dim.size();
		std::vector<int> id(ndim);
		for (int i=ndim-1; i>=0; --i){
			int ix = loc % dim[i];
			loc = (loc-ix)/dim[i];
			id[i]=ix;
		}
		return id;
	}

	/// A utility function for testing purposes. Fills the tensor with incremental integers. 
	void fill_sequence(){
		for(int i=0; i<vec.size(); ++i) vec[i]=i;
	}

	
	/// @brief generate a list of 1D indices corresponding to all points on the hyperplane 
	/// perpendicular to 'axis' located at index 'k' on the axis. The axis is specified
	/// as the index of the corresponding dimension, i.e., between [0, n-1].
	std::vector<int> plane(int axis, int k = 0){
		axis = dim.size()-1-axis;
		std::vector<int> locs;
		locs.reserve(nelem);
		for (int i=0; i<nelem; ++i){
			if (index(i)[axis] == 0) locs.push_back(i + k*offsets[axis]);
		}
		return locs;
	}

	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	void transform_dim(int loc, int axis, BinOp binary_op, std::vector<double> w){
		assert(w.size() == dim[dim.size()-1-axis]);
		
		axis = dim.size()-1-axis;
		int off = offsets[axis];
		
		for (int i=loc, count=0; count<dim[axis]; i+= off, ++count){
			vec[i] = binary_op(vec[i], w[count]);	// this order is important, because the operator may not be commutative
		}
		
	}

	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	void transform(int axis, BinOp binary_op, std::vector<double> w){
		std::vector<int> locs = plane(axis);
		for (int i=0; i<locs.size(); ++i){
			transform_dim(locs[i], axis, binary_op, w);
		}
	}
	
	
	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	double accumulate_dim(double v0, int loc, int axis, BinOp binary_op, std::vector<double> weights={}){
		assert(weights.size() == 0 || weights.size() == dim[dim.size()-1-axis]);
		
		axis = dim.size()-1-axis;
		int off = offsets[axis];
		
		double v = 0;
		for (int i=loc, count=0; count<dim[axis]; i+= off, ++count){
			double w = (weights.size()>0)? weights[count] : 1;
			v = binary_op(v, w*vec[i]);
		}
		
		return v;
	}


	// axis is counted from the right
	// [..., 2, 1, 0]
	//          ^
	//           axis
	template <class BinOp>
	Tensor<T> accumulate(T v0, int axis, BinOp binary_op, std::vector<double> weights={}){
		std::vector<int> dim_new = dim;
		dim_new.erase(dim_new.begin()+dim_new.size()-1-axis);
		Tensor<T> tens(dim_new);
		
		std::vector<int> locs = plane(axis);
		
		for (int i=0; i<locs.size(); ++i){
			tens.vec[i] = accumulate_dim(v0, locs[i], axis, binary_op, weights);
		}
		
		return tens;
	}


	Tensor<T> max_dim(int axis){
		T v0 = vec[1];
		return accumulate(v0, axis, [](T a, T b){return std::max(a,b);});
	}

	Tensor<T> avg_dim(int axis, std::vector<double> weights={}){
		Tensor tens = accumulate(0, axis, std::plus<T>(), weights);
		tens /= double(dim[dim.size()-1-axis]);
		return tens;
	}


	Tensor<T> repeat_inner(int n) const {
		std::vector<int> dim_new = dim;
		dim_new.push_back(n);
		
		Tensor<T> tout(dim_new);
		int count = 0;
		for (int i=0; i<nelem; ++i){
			for (int j=0; j<n; ++j){
				tout.vec[count++] = vec[i];
			}
		}

		return tout;
	}

	Tensor<T> repeat_outer(int n) const {
		std::vector<int> dim_new = dim;
		dim_new.insert(dim_new.begin(), n);
		
		Tensor<T> tout(dim_new);
		int count = 0;
		for (int j=0; j<n; ++j){
			for (int i=0; i<nelem; ++i){
				tout.vec[count++] = vec[i];
			}
		}

		return tout;
	}


	// operators
	public: 	
	// see https://stackoverflow.com/questions/4421706/what-are-the-basic-rules-and-idioms-for-operator-overloading/4421719#4421719
	template <class S>
	Tensor<T>& operator += (const Tensor<S>& rhs){
		assert(dim == rhs.dim);
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), std::plus<T>());
		return *this;
	}
	
	template <class S>
	Tensor<T>& operator -= (const Tensor<S>& rhs){
		assert(dim == rhs.dim);
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), std::minus<double>());
		return *this;
	}

	template <class S>
	Tensor<T>& operator *= (const Tensor<S>& rhs){
		assert(dim == rhs.dim);
		std::transform(vec.begin(), vec.end(), rhs.vec.begin(), vec.begin(), std::multiplies<T>());
		return *this;
	}

	template<class S>	
	Tensor<T>& operator += (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), [&s](const T& x){return x+s;});
		return *this;
	}

	template <class S>
	Tensor<T>& operator -= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), [&s](const T& x){return x-s;});
		return *this;
	}

	template <class S>
	Tensor<T>& operator *= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), [&s](const T& x){return x*s;});
		return *this;
	}

	template<class S>
	Tensor<T>& operator /= (S s){
		std::transform(vec.begin(), vec.end(), vec.begin(), [&s](const T& x){return x/s;});
		return *this;
	}

};

template<class T>
Tensor<T> operator + (Tensor<T> lhs, const Tensor<T>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs += rhs;
	return lhs;
}

template<class T>
Tensor<T> operator - (Tensor<T> lhs, const Tensor<T>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs -= rhs;
	return lhs;
}

template<class T>
Tensor<T> operator * (Tensor<T> lhs, const Tensor<T>& rhs){
	assert(lhs.dim == rhs.dim);
	lhs *= rhs;
	return lhs;
}

template<class T, class S>
Tensor<T> operator + (Tensor<T> lhs, S s){
	lhs += s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator - (Tensor<T> lhs, S s){
	lhs -= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator / (Tensor<T> lhs, S s){
	lhs /= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator * (Tensor<T> lhs, S s){
	lhs *= s;
	return lhs;
}

template<class T, class S>
Tensor<T> operator + (S s, Tensor<T> t){	// TODO: can passing be ref be used for these?
	return t+s;
}

template<class T, class S>
Tensor<T> operator - (S s, Tensor<T> t){
	return t-s;
}

template<class T, class S>
Tensor<T> operator * (S s, Tensor<T> t){
	return t*s;
}




#endif
