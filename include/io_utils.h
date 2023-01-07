#ifndef PSPM_IO_UTILS_H_
#define PSPM_IO_UTILS_H_

#include <ostream>
#include <istream>
#include <tuple>
#include <vector>

// print a vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	os << '\n';
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


// print a tuple (requires C++17)
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...> t) {
    apply([&](auto&&... args) {
        ((os << args << " "), ...);
    }, t);
    return os;
}



#endif

