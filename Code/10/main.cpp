#include "nr3.h"
#include "utils/ode.h"
#include <format>
#include <print>

VecDoub derivs(const Doub x, VecDoub_I &y) {
  (void)x; // Unused
  VecDoub_O dydx(2);
  dydx[0] = y[0] * y[1];
  dydx[1] = -(y[0] * y[0]);
  return dydx;
}

int main() {

  {
    std::println("\n1nd order (Euler method):");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    print_ode_table(low, high, 5, y, derivs, ODEMethod::Euler, 1e-5);
    // for (int n = 5; n <= 1000; n *= 2) {
    //   auto y_n = y;
    //   y_n = first_order_runge_kuttea_method(low, high, n, y_n, derivs);
    //   util::print(y_n, std::format("1st order with h={}", (high - low) / n));
    // }
  }

  {
    std::println("\n2nd order (Midpoint method):");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    print_ode_table(low, high, 5, y, derivs, ODEMethod::Midpoint, 1e-5);
    // for (int n = 5; n <= 1000; n *= 2) {
    //   const auto res =
    //       second_order_runge_kuttea_method(low, high, n, y, derivs);
    //   util::print(res,
    //               std::format("2st order with h={}", (high - low) /
    //               (double)n));
    // }
  }

  {
    std::println("\n4nd order:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    print_ode_table(low, high, 5, y, derivs, ODEMethod::FORTH_ORDER, 1e-5);
    // for (int n = 5; n <= 1000; n *= 2) {
    //   const auto res =
    //       fourth_order_runge_kuttea_method(low, high, n, y, derivs);
    //   util::print(res,
    //               std::format("4st order with h={}", (high - low) /
    //               (double)n));
    // }
  }

  {
    std::println("\nTrapezoidal method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    print_ode_table(low, high, 5, y, derivs, ODEMethod::TRAPEZOIDAL, 1e-5);
    // for (int n = 5; n <= 1000; n *= 2) {
    //   auto y_n = y;
    //   VecDoub res = trapezoidal(low, high, n, y_n, derivs);
    //   util::print(
    //       res, std::format("Trapezoidal with h={}", (high - low) /
    //       (double)n));
    // }
  }

  {
    // INFO: Very unstable method...
    std::println("\nLeap Frog method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    print_ode_table(low, high, 5, y, derivs, ODEMethod::LEAP_FROG, 1e-5);
    // for (int n = 5; n <= 1000; n *= 2) {
    //   auto y_n = y;
    //   VecDoub res = leap_frog_method(low, high, n, y_n, derivs);
    //   util::print(res,
    //               std::format("Leap Frog with h={}", (high - low) /
    //               (double)n));
    // }
  }
}
