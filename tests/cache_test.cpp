#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
using namespace std;

#include "iterator_set.h"

double f(double x){
	return (x);
}


int main(){


	int n = 10000;
	vector <double> state(5*n);
	
	{
	// x x x x y y y y a b c a b c a b c a b c 
	auto start = chrono::steady_clock::now();
	for (int k=0; k<1000; ++k){
		auto ix = state.begin();
		auto iy = next(state.begin(), n);
		auto ia = next(state.begin(), 2*n + 0);
		auto ib = next(state.begin(), 2*n + 1);
		auto ic = next(state.begin(), 2*n + 2);

		for (int i=0; i<n; ++i){
			double x = double(i)/(n-1); 
			*ix = x;
			*iy = f(-x);

			*ia = f(-x);
			*ib = f(-x)*10;
			*ic = f(-x)*20;

			++ix; ++iy; 
			ia += 3;
			ib += 3;
			ic += 3;
		}
	}
	auto end = chrono::steady_clock::now();
	cout << "t xx yy abc abc  = " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	}

	{
	// x x x x y y y y a b c a b c a b c a b c 
	auto start = chrono::steady_clock::now();
	for (int k=0; k<1000; ++k){
		IteratorSet<vector<double>::iterator> iset(state.begin(), {"x","y","a","b","c"}, n, {n,n,n,n,n}, {1,1,1,1,1});
		auto& iters = iset.get();
		auto ix = iters[0];
		auto iy = iters[1];
		auto ia = iters[2];
		auto ib = iters[3];
		auto ic = iters[4];

		int i=0;
		for (; i<n; ++i){
			double x = double(i)/(n-1); 
			*ix = x;
			*iy = f(-x);

			*ia = f(-x);
			*ib = f(-x)*10;
			*ic = f(-x)*20;
		
			//++iset;	
			//++ix; ++iy; ++ia;++ib;++ic;
			ix+=iset.strides[0];
			iy+=iset.strides[1];
			ia+=iset.strides[2];
			ib+=iset.strides[3];
			ic+=iset.strides[4];
		}
	}
	auto end = chrono::steady_clock::now();
	cout << "T xx yy abc abc  = " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	//iset.printHeader();
	//for (iset.begin(); !iset.end(); ++iset) iset.print();
	}


	{
	// x x x x y y y y a a a a b b b b c c c c 
	auto start = chrono::steady_clock::now();
	for (int k=0; k<1000; ++k){
		auto ix = state.begin();
		auto iy = next(state.begin(), n);
		auto ia = next(state.begin(), 2*n);
		auto ib = next(state.begin(), 3*n);
		auto ic = next(state.begin(), 4*n);

		for (int i=0; i<n; ++i){
			double x = double(i)/(n-1);
 
			*ix = x;
			*iy = f(-x);

			*ia = f(-x);
			*ib = f(-x)*10;
			*ic = f(-x)*20;

			++ix; ++iy; 
			++ia; //advance(ia, 1);
			++ib; //advance(ib, 1);
			++ic; //advance(ic, 1);
		}
	}
	auto end = chrono::steady_clock::now();
	cout << "t xx yy aa bb cc = " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	IteratorSet<vector<double>::iterator> iset(state.begin(), {"x","y","a","b","c"}, n, {n,n,n,n,n}, {1,1,1,1,1});
	iset.print(10);
	}

	{
	// x x x x y y y y a a a a b b b b c c c c 
	auto start = chrono::steady_clock::now();
	for (int k=0; k<1000; ++k){
		IteratorSet<vector<double>::iterator> iset(state.begin(), {"x","y","a","b","c"}, n, {n,n,n,n,n}, {1,1,1,1,1});
		auto& iters = iset.get();
		auto ix = iters[0];
		auto iy = iters[1];
		auto ia = iters[2];
		auto ib = iters[3];
		auto ic = iters[4];

		int i=0;
		for (; i<n; ++i){
			double x = double(i)/(n-1); 
			*ix = x;
			*iy = f(-x);

			*ia = f(-x);
			*ib = f(-x)*10;
			*ic = f(-x)*20;
			
			//++ix; ++iy; ++ia;++ib;++ic;
			ix+=iset.strides[0];
			iy+=iset.strides[1];
			ia+=iset.strides[2];
			ib+=iset.strides[3];
			ic+=iset.strides[4];

		}
	}
	auto end = chrono::steady_clock::now();
	cout << "T xx yy aa bb cc = " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	IteratorSet<vector<double>::iterator> iset(state.begin(), {"x","y","a","b","c"}, n, {n,n,n,n,n}, {1,1,1,1,1});
	iset.print(10);	
	}
	
	{
	// x y a b c x y a b c...  
	auto start = chrono::steady_clock::now();
	for (int k=0; k<1000; ++k){
		auto ix = state.begin();

		for (int i=0; i<n; ++i){
			double x = double(i)/(n-1);
			*ix++ = x;
			*ix++ = f(-x);
			*ix++ = f(-x);
			*ix++ = f(-x)*10;
			*ix++ = f(-x)*20;
			
		}
	}
	auto end = chrono::steady_clock::now();
	cout << "t xyabc xyabc    = " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
	}


	//// print
	//auto ix = state.begin();
	//auto iy = next(state.begin(), n);
	//auto ia = next(state.begin(), 2*n + 0);
	//auto ib = next(state.begin(), 3*n + 0);
	//auto ic = next(state.begin(), 4*n + 0);
	//for (int i=0; i<n; ++i){
	//    cout << *ix << "\t" << *iy << "\t" << *ia << "\t" << *ib << "\t" << *ic << "\n";
	//    ++ix; ++iy; 
	//    advance(ia, 1);
	//    advance(ib, 1);
	//    advance(ic, 1);
	//}	
	//cout << endl;



}
