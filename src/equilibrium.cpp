// This equilibrium solver is based on the python prototype implemented by Aseem Godbole on 1 Nov 2025

#include <cmath>
#include <iostream>
#include "elementwise.h"
#include "index_utils.h"
#include "solver.h"
using namespace std;


// Environment calculation function - replace suitably with Solver::updateEnv()
// Elementwise<double> E_func(const Elementwise<double>& u, const Elementwise<double>& x, const Elementwise<double>& z) {
//     Elementwise<double> Ez(z.size(), 0.0);
    
//     for (size_t j = 0; j < z.size(); ++j) {
//         Elementwise<double> x_filtered, u_filtered;
//         for (size_t i = 0; i < x.size(); ++i) {
//             if (x[i] >= z[j]) {
//                 x_filtered.push_back(x[i]);
//                 u_filtered.push_back(u[i] * Ac(x[i]));
//             }
//         }
//         if (!x_filtered.rawdata().empty()) {
//             Ez[j] = trapz(u_filtered, x_filtered);
//         }
//     }
    
//     return Ez * E0;  // Now uses operator* directly
// }

// Elementwise<double> int_1_over_g(const Elementwise<double>& x, const Elementwise<double>& Ez) {
//     return cumtrapz(x.transform([](double xi) { return 1.0/g(xi); }), x);
// }

// Implement this function
// std::pair<Elementwise<double>, Elementwise<double>> Solver::solve_U_equil() {
//     Elementwise<double> x = log_seq(1.0, 1e6, 50);
//     Elementwise<double> z = log_seq(1.0, 1e6, 200);
    
//     // Initialize u = x^(-phi_g)
//     Elementwise<double> u = x.transform([](double xi) { return pow(xi, -phi_g); });
//     double norm = trapz(u, x);
//     u = u / norm;  // Now uses operator/ directly
    
//     Elementwise<double> Ez = E_func(u, x, z);
//     int cnt = 0;
    
//     for (int it = 0; it < 2000; ++it) {
//         Ez = E_func(u, x, z);
//         cnt++;
        
//         Elementwise<double> gx = x.transform(g);

//         Elementwise<double> I = int_1_over_g(x, Ez);
        
//         // Safe clipping of u
//         Elementwise<double> u_safe = u.transform([](double ui) { 
//             return std::clamp(ui, 1e-20, 1e3); 
//         });
        
//         // Calculate Bx and B
//         Elementwise<double> Bx = x.transform([&Ez](double xi) { return betax(xi, Ez); });
//         double B = trapz(Bx * u_safe, x);  // Now uses operator* directly
        
//         // Calculate new u
//         Elementwise<double> u_new = I.transform([](double Ii) { return exp(-myu0 * Ii); }) / gx * B;
        
//         // Check convergence
//         double max_diff = 0.0;
//         for (size_t i = 0; i < x.size(); ++i) {
//             max_diff = std::max(max_diff, std::abs(u_new[i] - u[i]));
//         }
        
//         if (max_diff < 1e-6) {
//             std::cout << "Converged after " << cnt << " iterations\n";
//             return {u_new, Ez};
//         }
//         u = u_new;
//     }
//     return {u, Ez};
// }

// int main() {
//     auto [u_final, Ez_final] = solve_u();
//     std::cout << u_final << std::endl;
//     return 0;
// }
