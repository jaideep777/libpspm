#pragma once
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <iostream>

template<typename T>
class Elementwise {
private:
    std::vector<T> data;

public:
    // Constructors
    Elementwise() = default;
    explicit Elementwise(size_t size) : data(size) {}
    Elementwise(size_t size, const T& value) : data(size, value) {}
    Elementwise(std::initializer_list<T> init) : data(init) {}
    
    // Convert from std::vector
    explicit Elementwise(const std::vector<T>& v) : data(v) {}
    
    // Copy/conversion constructor
    template<typename U>
    Elementwise(const Elementwise<U>& other) {
        data.resize(other.size());
        std::copy(other.begin(), other.end(), begin());
    }

    // Assignment
    Elementwise& operator=(const Elementwise& other) {
        if (this != &other) data = other.data;
        return *this;
    }

    // Vector access
    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }
    
    // Size operations
    size_t size() const { return data.size(); }
    void resize(size_t n) { data.resize(n); }
    void push_back(const T& value) { data.push_back(value); }
    
    // Iterators
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }

    // Access to data
    const std::vector<T>& rawdata() const { return data; }

    // Compound assignment operators
    template<typename U>
    Elementwise& operator+=(const Elementwise<U>& rhs) {
        std::transform(begin(), end(), rhs.begin(), begin(), std::plus<>());
        return *this;
    }

    template<typename U>
    Elementwise& operator-=(const Elementwise<U>& rhs) {
        std::transform(begin(), end(), rhs.begin(), begin(), std::minus<>());
        return *this;
    }

    template<typename U>
    Elementwise& operator*=(const Elementwise<U>& rhs) {
        std::transform(begin(), end(), rhs.begin(), begin(), std::multiplies<>());
        return *this;
    }

    template<typename U>
    Elementwise& operator/=(const Elementwise<U>& rhs) {
        std::transform(begin(), end(), rhs.begin(), begin(), std::divides<>());
        return *this;
    }
 
    template<typename U>
    Elementwise& operator*=(U scalar) {
        std::transform(begin(), end(), begin(),
                      [scalar](const T& x) { return x * scalar; });
        return *this;
    }

    template<typename U>
    Elementwise& operator/=(U scalar) {
        std::transform(begin(), end(), begin(),
                      [scalar](const T& x) { return x / scalar; });
        return *this;
    }

    template<typename U>
    Elementwise& operator+=(U scalar) {
        std::transform(begin(), end(), begin(),
                      [scalar](const T& x) { return x + scalar; });
        return *this;
    }

    template<typename U>
    Elementwise& operator-=(U scalar) {
        std::transform(begin(), end(), begin(),
                      [scalar](const T& x) { return x - scalar; });
        return *this;
    }

    // Transform
    template<typename Func>
    Elementwise<T> transform(Func f) const {
        Elementwise<T> result(size());
        std::transform(begin(), end(), result.begin(), f);
        return result;
    }

    // Magnitude (L2 norm)
    T magnitude() const {
        return std::sqrt(std::inner_product(begin(), end(), begin(), T()));
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const Elementwise<T>& v) {
        os << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) os << ", ";
            os << v[i];
        }
        os << "]";
        return os;
    }
};

// Elementwise vector multiplication
template<typename T, typename U>
Elementwise<T> operator*(Elementwise<T> lhs, const Elementwise<U>& rhs) {
    lhs *= rhs;
    return lhs;
}

// Add elementwise division operator after other operators
template<typename T, typename U>
Elementwise<T> operator/(Elementwise<T> lhs, const Elementwise<U>& rhs) {
    lhs /= rhs;
    return lhs;
}

// Standalone operators
template<typename T, typename U>
Elementwise<T> operator+(Elementwise<T> lhs, const Elementwise<U>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T, typename U>
Elementwise<T> operator-(Elementwise<T> lhs, const Elementwise<U>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, typename U>
Elementwise<T> operator*(Elementwise<T> v, U scalar) {
    v *= scalar;
    return v;
}

template<typename T, typename U>
Elementwise<T> operator/(Elementwise<T> v, U scalar) {
    v /= scalar;
    return v;
}

template<typename T, typename U>
Elementwise<T> operator+(Elementwise<T> v, U scalar) {
    v += scalar;
    return v;
}

template<typename T, typename U>
Elementwise<T> operator-(Elementwise<T> v, U scalar) {
    v -= scalar;
    return v;
}


// Utility functions
template<typename T>
Elementwise<T> normalize(const Elementwise<T>& v, T tol = T()) {
    T mag = v.magnitude() + tol;
    return v / mag;
}

template<typename T>
T dot(const Elementwise<T>& lhs, const Elementwise<T>& rhs) {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), T());
}

template<typename T>
Elementwise<T> reciprocal(const Elementwise<T>& v) {
    Elementwise<T> res(v.size());
    std::transform(v.begin(), v.end(), res.begin(),
                   [](const T& x) { return 1.0 / x; });
    return res;
}
