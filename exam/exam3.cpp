#include "nr3.h"
#include "utils/ode.h"
#include "utils/utilities.h"
#include <algorithm>
#include <math.h>
#include <print>

double a_max = 4.0;
double v_des = 25.0;
double D_0 = 50.0;
double T_react = 1.5;
double a_com = 2.0;

double X_F(double t) { return 250 + 15 * t - 5 * sqrt(1 + pow(t, 2)); }
double X_F_prime(double t) { return 15 - (5 * t) / sqrt(1 + pow(t, 2)); }

VecDoub derivs(const Doub t, VecDoub_I &x) {
  double x1 = x[0];
  double x2 = x[1];
  VecDoub_O dxdt(2);
  dxdt[0] = x2;
  dxdt[1] =
      a_max *
      (1 - pow((x2 / v_des), 4) -
       pow((D_0 + std::max(0.0, x2 * T_react +
                                    (x2 * (x2 - X_F_prime(t))) / (2 * a_com))) /
               (X_F(t) - x1),
           2));
  return dxdt;
}

// Midpoint method
VecDoub second_order_runge_kuttea_method_plot(double low, double high,
                                              int steps, VecDoub_I &x,
                                              VecDoub derivs(const Doub t,
                                                             VecDoub_I &x)) {
  std::vector<std::vector<double>> plotting;
  const double h = (high - low) / (double)steps;
  VecDoub x_n = x;
  plotting.push_back({low, x_n[0], x_n[1], (X_F(low) - x_n[0])});
  for (double t_n = low; t_n < high; t_n += h) {
    f_comps_current += 2;
    auto k1 = h * derivs(t_n, x_n);
    auto k2 = h * derivs(t_n + 0.5 * h, x_n + 0.5 * k1);
    auto x_n_next = x_n + k2;
    x_n = x_n_next;
    plotting.push_back({t_n + h, x_n[0], x_n[1], (X_F(t_n + h) - x_n[0])});
  }

  std::ofstream file;
  file.open("exam3_results.txt");
  for (const auto &vec : plotting) {
    std::println(file, "{}\t{}\t{}\t{}", vec[0], vec[1], vec[2], vec[3]);
  }
  file.close();

  return x_n;
}

int main() {
  {
    // 2)
    std::println("\n{:/^100}", " i ");

    double initial_vals[2] = {0.0, 15.0};
    VecDoub x(2, initial_vals);
    auto res = derivs(-10, x);
    util::print(res, "res");
  }

  {
    // 3)
    std::println("\n{:/^100}", " ii ");

    double initial_vals[2] = {0.0, 15.0};
    VecDoub x(2, initial_vals);
    auto res = second_order_runge_kuttea_method_plot(-10, 10, 80, x, derivs);
    util::print(res, "res");
  }

  {
    // 4)
    std::println("\n{:/^100}", " iii ");
  }

  {
    // 5)
    std::println("\n{:/^100}", " iv ");
    double initial_vals[2] = {0.0, 15.0};
    VecDoub x(2, initial_vals);
    print_ode_table(-10, 10, 20, x, derivs, ODEMethod::Midpoint, 2e-5);
  }
}
