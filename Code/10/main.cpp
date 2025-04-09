#include "nr3.h"
#include "utilities.h"
#include <format>
#include <print>

// Euler method
// Midpoint method
// Trapezoidal method
// 4th order Runge-Kutta method
// Leap-Frog method
// for
// 0 ≤ x ≤ 10,  n = 5, 10, 20, 40 . . .

VecDoub derivs(const Doub x, VecDoub_I &y) {
  (void)x; // Unused
  VecDoub_O dydx(2);
  dydx[0] = y[0] * y[1];
  dydx[1] = -(y[0] * y[0]);
  return dydx;
}

VecDoub first_order_runge_kuttea_method(double x, VecDoub_I &y, const Doub h,
                                        VecDoub derivs(const Doub x,
                                                       VecDoub_I &y)) {
  // Euler method
  return y + h * derivs(x, y);
}

VecDoub second_order_runge_kuttea_method(double x, VecDoub_I &y, const Doub h,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  // Midpoint method
  auto k1 = h * derivs(x, y);
  auto k2 = h * derivs(x + 0.5 * h, y + 0.5 * k1);
  return y + k2;
}

VecDoub third_order_runge_kuttea_method(double x0, VecDoub_I &y, const Doub h,
                                        VecDoub derivs(const Doub x,
                                                       VecDoub_I &y)) {
  // Trapezoidal method
}

VecDoub fourth_order_runge_kuttea_method(double x, VecDoub_I &y, const Doub h,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {}

VecDoub leap_frog_method(double x, VecDoub_I &y, const Doub h,
                         VecDoub derivs(const Doub x, VecDoub_I &y)) {
  return y + 2 * h * derivs(x, y);
}

int main() {
  {
    std::println("\n1nd order:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 1000; n *= 2) {
      double h = (high - low) / n;
      auto y_n = y;
      for (int i = 0; i < n; i++) {
        VecDoub res = first_order_runge_kuttea_method(0, y_n, h, derivs);
        result.push_back(res);
        y_n = res;
      }
      util::print(y_n, std::format("1st order with h={}", h));
    }
  }

  {
    std::println("\n2nd order:");
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 100; n *= 2) {
      double h = (high - low) / n;
      auto y_n = y;
      for (int i = 0; i < n; i++) {
        VecDoub res = second_order_runge_kuttea_method(0, y_n, h, derivs);
        result.push_back(res);
        y_n = res;
      }
      util::print(y_n, std::format("1st order with h={}", h));
    }
  }
  //
  // {
  //   double vals[2] = {1.0, 1.0};
  //   const VecDoub y(2, vals);
  //   double low = 0.0;
  //   double high = 10.0;
  //   std::vector<VecDoub> result;
  //   for (int n = 5; n <= 100; n *= 2) {
  //     double h = (high - low) / n;
  //     auto y_n = y;
  //     for (int i = 0; i < n; i++) {
  //       VecDoub res = third_order_runge_kuttea_method(0, y_n, h, derivs);
  //       result.push_back(res);
  //       y_n = res;
  //     }
  //     util::print(y_n, std::format("1st order with h={}", h));
  //   }
  // }
  //
  // {
  //   double vals[2] = {1.0, 1.0};
  //   const VecDoub y(2, vals);
  //   double low = 0.0;
  //   double high = 10.0;
  //   std::vector<VecDoub> result;
  //   for (int n = 5; n <= 100; n *= 2) {
  //     double h = (high - low) / n;
  //     auto y_n = y;
  //     for (int i = 0; i < n; i++) {
  //       VecDoub res = fourth_order_runge_kuttea_method(0, y_n, h, derivs);
  //       result.push_back(res);
  //       y_n = res;
  //     }
  //     util::print(y_n, std::format("1st order with h={}", h));
  //   }
  // }
  {
    double vals[2] = {1.0, 1.0};
    const VecDoub y(2, vals);
    double low = 0.0;
    double high = 10.0;
    std::vector<VecDoub> result;
    for (int n = 5; n <= 1000; n *= 2) {
      double h = (high - low) / n;
      auto y_n = y;
      for (int i = 0; i < n; i++) {
        if (i == 0) {
          VecDoub res = first_order_runge_kuttea_method(0, y_n, h, derivs);
          result.push_back(res);
          y_n = res;
        } else {
          VecDoub res = leap_frog_method(0, y_n, h, derivs);
          result.push_back(res);
          y_n = res;
        }
      }
      util::print(y_n, std::format("leap_frog with h={}", h));
    }
  }
}
