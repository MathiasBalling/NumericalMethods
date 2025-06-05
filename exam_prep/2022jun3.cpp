#include "nr3.h"
#include "utils/ode.h"
#include "utils/utilities.h"
#include <print>
const double a1 = 1.0, b1 = 1.0, a2 = 0.25, b2 = 0.25;
const double Y1 = 0.0, Y2 = 0.0, g1 = 2.0, g2 = 1.0;

double f1(double x) { return x; }
double f2(double x) { return 3 * x * (1 - pow(x, 2)); }

VecDoub derivs(const Doub t, VecDoub_I &v) {
  const double z1 = v[0]; // y_1(t)
  const double z2 = v[1]; // y_1'(t)
  const double z3 = v[2]; // y_2(t)
  const double z4 = v[3]; // y_2'(t)
  const double z5 = v[4]; // x(t)

  VecDoub_O dzdt(5);
  dzdt[0] = z2;
  dzdt[1] = a1 * (b1 * (g1 - z1) - z2) + f1(z5);
  dzdt[2] = z4;
  dzdt[3] = a2 * (b2 * (g2 - z3) - z4) + f2(z5);
  dzdt[4] = -z5;
  return dzdt;
}

// Midpoint method
void second_order_runge_kuttea_method_plot(double low, double high, int steps,
                                           VecDoub_I &y,
                                           VecDoub derivs(const Doub x,
                                                          VecDoub_I &y)) {
  std::vector<VecDoub> y_n_trajectory;
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  y_n_trajectory.push_back(y_n);
  for (double x_n = low; x_n < high; x_n += h) {
    f_comps_current += 2;
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto y_n_next = y_n + k2;
    y_n = y_n_next;
    y_n_trajectory.push_back(y_n);
  }
  std::ofstream file;
  file.open("2022jun3_results.txt");
  for (const auto &vec : y_n_trajectory) {
    std::println(file, "{}\t{}", vec[0], vec[2]);
  }
  file.close();
}

int main() {
  {
    // 1)
    // $$
    // \begin{aligned}
    // y_1''(t) &= a_1\left\{b_1[g_1 - y_1(t)] - y_1'(t)\right\} + f_1(x(t)),
    // \quad y_1(0) = Y_1, \quad y_1'(0) = 0 \\
    // y_2''(t) &= a_2\left\{b_2[g_2 - y_2(t)] - y_2'(t)\right\} + f_2(x(t)),
    // \quad y_2(0) = Y_2, \quad y_2'(0) = 0 \\ x'(t) &= -x(t), \quad x(0) = 1
    // \end{aligned}
    // $$
    //
    // ### Step 1: Define new variables to reduce to first-order
    //
    // * $z_1(t) = y_1(t)$
    // * $z_2(t) = y_1'(t)$
    // * $z_3(t) = y_2(t)$
    // * $z_4(t) = y_2'(t)$
    // * $z_5(t) = x(t)$
    //
    // ### Step 2: Express as first-order system
    //
    // $$
    // \begin{aligned}
    // z_1'(t) &= z_2(t) \\
    // z_2'(t) &= a_1\left[b_1(g_1 - z_1(t)) - z_2(t)\right] + f_1(z_5(t)) \\
    // z_3'(t) &= z_4(t) \\
    // z_4'(t) &= a_2\left[b_2(g_2 - z_3(t)) - z_4(t)\right] + f_2(z_5(t)) \\
    // z_5'(t) &= -z_5(t)
    // \end{aligned}
    // $$
    // ### Step 3: Initial Conditions
    //
    // $$
    // \begin{aligned}
    // z_1(0) &= Y_1 \\
    // z_2(0) &= 0 \\
    // z_3(0) &= Y_2 \\
    // z_4(0) &= 0 \\
    // z_5(0) &= 1
    // \end{aligned}
    // $$
  }

  {
    // 2)
    // std::println("\n{:/^80}", " ii ");
    double initial_vals[5] = {Y1, 0, Y2, 0.0, 1.0};
    const VecDoub y(5, initial_vals);
    double t_low = 0.0;
    double t_high = 20.0;
    double h = 0.001;
    int steps = (t_high - t_low) / h;

    second_order_runge_kuttea_method_plot(t_low, t_high, steps, y, derivs);
  }

  {
    // 3)
    std::println("\n{:/^80}", " iii ");
    double initial_vals[5] = {Y1, 0, Y2, 0.0, 1.0};
    const VecDoub y(5, initial_vals);
    double t_low = 0.0;
    double t_high = 5.0;

    print_ode_table(t_low, t_high, 5, y, derivs, ODEMethod::Midpoint, 1e-4);
  }
}
