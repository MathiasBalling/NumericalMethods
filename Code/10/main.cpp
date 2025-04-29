#include "multi_roots_table.h"
#include "nr3.h"
#include "roots.h"
#include "roots_multidim.h"
#include "utilities.h"
#include <format>
#include <numbers>
#include <print>

VecDoub derivs(const Doub x, VecDoub_I &y) {
  (void)x; // Unused
  VecDoub_O dydx(2);
  dydx[0] = y[0] * y[1];
  dydx[1] = -(y[0] * y[0]);
  return dydx;
}

VecDoub
first_order_runge_kuttea_method(double low, double high, int steps, VecDoub_I y,
                                VecDoub derivs(const Doub x, VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;

  VecDoub y_n = y;
  for (double x = low; x < high; x += h) {
    auto y_n_next = y_n + h * derivs(x, y_n);
    y_n = y_n_next;
  }
  // Euler method
  return y_n;
}

VecDoub second_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &y,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  // Midpoint method
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto y_n_next = y_n + k2;
    y_n = y_n_next;
  }
  return y_n;
}

VecDoub third_order_runge_kuttea_method(double low, double high, int steps,
                                        VecDoub_I &y,
                                        VecDoub derivs(const Doub x,
                                                       VecDoub_I &y)) {
  // Trapezoidal method
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    // Euler step
    auto y_n_next = y_n + h * derivs(x_n, y_n);

    auto func = [=](VecDoub y_guess) {
      return y_guess - y_n -
             h / 2 * (derivs(x_n, y_n) + derivs(x_n + h, y_guess));
    };

    bool check;
    newton_multi(y_n_next, check, func, 200, 1.0e-7);
    y_n = y_n_next;
  }

  return y_n;
}

VecDoub leap_frog_method(double low, double high, int steps, VecDoub_I &y,
                         VecDoub derivs(const Doub x, VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub y_n_last = y + h * derivs(low, y);
  VecDoub y_n = y_n_last + h * derivs(low + h, y_n_last);
  for (double x_n = low + h + h; x_n < high; x_n += h) {
    VecDoub y_n_next = y_n_last + 2 * h * derivs(x_n, y_n);
    y_n_last = y_n;
    y_n = y_n_next;
  }
  return y_n;
}

int main() {
  {
    std::println("\n1nd order Euler method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 20000; n *= 2) {
      auto y_n = y;
      y_n = first_order_runge_kuttea_method(low, high, n, y_n, derivs);
      util::print(y_n, std::format("1st order with h={}", (high - low) / n));
    }
  }

  {
    std::println("\n2nd order Midpoint method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 1000; n *= 2) {
      const auto res =
          second_order_runge_kuttea_method(low, high, n, y, derivs);
      util::print(res,
                  std::format("2st order with h={}", (high - low) / (double)n));
    }
  }

  {
    std::println("\n3nd order Trapezoidal method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 1000; n *= 2) {
      auto y_n = y;
      VecDoub res = third_order_runge_kuttea_method(low, high, n, y_n, derivs);
      util::print(res,
                  std::format("3st order with h={}", (high - low) / (double)n));
    }
  }

  {
    std::println("\n4nd order Leap Frog method:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 1000; n *= 2) {
      auto y_n = y;
      VecDoub res = leap_frog_method(low, high, n, y_n, derivs);
      util::print(res,
                  std::format("4st order with h={}", (high - low) / (double)n));
    }
  }
}
