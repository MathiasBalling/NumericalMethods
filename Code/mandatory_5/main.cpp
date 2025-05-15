#include "utils/bvp.h"
#include <print>
// Problem:
// V_mm(u)  = (48(V(u)^3+2u^3*V'(u))(2u^2-V(u)^2*V'(u)^2))/(1+64u^6+16V(u)^6
// V(-0.9)=-0.85
// V(0.8)=-0.9

// V'' = F
double F(double v_prime, double v, double u) {
  return 48 * (pow(v, 3) + 2 * pow(u, 3) * v_prime) *
         (2 * pow(u, 2) - pow(v, 2) * pow(v_prime, 2)) /
         (1 + 64 * pow(u, 6) + 16 * pow(v, 6));
};

// d/dV F = F_v
double F_v(double v_prime, double v, double u) {
  return (48 *
              (3 * pow(v, 2) * (2 * pow(u, 2) - pow(v, 2) * pow(v_prime, 2)) +
               (pow(v, 3) + 2 * pow(u, 3) * v_prime) *
                   (-2 * v * pow(v_prime, 2))) *
              (1 + 64 * pow(u, 6) + 16 * pow(v, 6)) -
          48 * (pow(v, 3) + 2 * pow(u, 3) * v_prime) *
              (2 * pow(u, 2) - pow(v, 2) * pow(v_prime, 2)) *
              (96 * pow(v, 5))) /
         pow(1 + 64 * pow(u, 6) + 16 * pow(v, 6), 2);
};

// d/d(V') F = F_v_prime
double F_v_prime(double v_prime, double v, double u) {
  return 48 *
         (2 * pow(u, 3) * (2 * pow(u, 2) - pow(v, 2) * pow(v_prime, 2)) +
          (pow(v, 3) + 2 * pow(u, 3) * v_prime) * (-2 * pow(v, 2) * v_prime)) /
         (1 + 64 * pow(u, 6) + 16 * pow(v, 6));
};

int main() {
  {
    double a = -0.9;
    double b = 0.8;
    double alpha = -0.85;
    double beta = -0.9;
    int N = 2;
    // u((-0.9+0.8)/2) = ?
    double target_u = (a + b) / 2.0;
    finite_difference_method_table(N, target_u, a, b, alpha, beta, F, F_v,
                                   F_v_prime, 1e-6, 10000);
  }
  {
    double a = -0.9;
    double b = 0.8;
    double alpha = -0.85;
    double beta = -0.9;
    int N = 2;
    // u((-0.9+0.8)/2) = ?
    double target_u = (a + b) / 2.0;
    const auto res = finite_difference_method(N, target_u, a, b, alpha, beta, F,
                                              F_v, F_v_prime, 1e-6, 10000);
    const auto final_u = res.back();
    // Write the results to a file for plotting
    std::ofstream file;
    file.open("results.txt");
    const auto N_final = final_u.size();
    const auto h_final = (b - a) / (double)N_final;
    for (int i = 0; i < N_final; i++) {
      std::println(file, "{}\t{}", a + i * h_final, final_u[i]);
    }
    file.close();
  }
}
