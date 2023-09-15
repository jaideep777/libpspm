#ifndef VECTOR_IO_UTILS_H_
#define VECTOR_IO_UTILS_H_

#include <ostream>
#include <istream>
#include <tuple>
#include <vector>
#include <array>
#include <iomanip>

// print a vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	//os << std::setprecision(12);
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	os << '\n'; // FIXME: remove this newline and insert newline in every save() function
	return os;
}


// read a vector via ofstream
// reads size, "|", discards "|",
// then resizes the vector and reads values
template <class T>
std::istream& operator >> (std::istream &is, std::vector<T> &v) {
	int n; 
	std::string s;
	is >> n >> s; // s contains "|"
	v.resize(n);
	for (int i=0; i<n; ++i) is >> v[i];
	return is;
}


// print an array via ofstream
// prints: size | v1 v2 v3 ...
template <class T, size_t n>
std::ostream& operator << (std::ostream &os, const std::array<T,n> &v) {
	//os << std::setprecision(12);
	for (const auto &x : v) {
		os << x << ' ';
	}
	// os << '\n'; // FIXME: remove this newline and insert newline in every save() function
	return os;
}


// read an array via ofstream
// reads size, "|", discards "|",
// then resizes the vector and reads values
template <class T, size_t n>
std::istream& operator >> (std::istream &is, std::array<T,n> &v) {
	for (int i=0; i<n; ++i) is >> v[i];
	return is;
}


// print a tuple (requires C++17)
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...> t) {
    apply([&](auto&&... args) {
        ((os << args << " "), ...);
    }, t);
    return os;
}



#endif

