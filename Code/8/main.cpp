#include <print>
// #include "quadrature.h"

template <class T>
auto midpoint(T &func, double limit_low, double limit_high, int its = 100)
    -> double {
  // TODO: Add table printing
  double step = (limit_high - limit_low) / its;
  double result = 0.0;
  for (int i = 0; i < its; ++i) {
    double x1 = limit_low + i * step;
    double x2 = x1 + step;
    result += func((x1 + x2) / 2) * step;
  }
  return result;
}

template <class T>
auto trapezoidal(T &func, double limit_low, double limit_high, int its = 100)
    -> double {
  // TODO: Add table printing
  double step = (limit_high - limit_low) / its;
  double result = 0.5 * (func(limit_low) + func(limit_high));
  for (int i = 1; i < its; ++i) {
    double x = limit_low + i * step;
    result += func(x);
  }
  result *= step;
  return result;
}

template <class T>
auto simpson(T &func, double limit_low, double limit_high, int its = 100)
    -> double {
  // TODO: Add table printing
  if (its % 2 != 0) {
    ++its; // Ensure even number of intervals
  }
  double step = (limit_high - limit_low) / its;
  double result = func(limit_low) + func(limit_high);
  for (int i = 1; i < its; ++i) {
    double x = limit_low + i * step;
    result += func(x) * (i % 2 == 0 ? 2 : 4);
  }
  result *= step / 3.0;
  return result;
}

//
auto eqn_1(double x) -> double { return cos(x * x) * exp(-x); }
auto eqn_2(double x) -> double { return sqrt(x) * cos(x * x) * exp(-x); }
auto eqn_3(double x) -> double { return (1 / sqrt(x)) * cos(x * x) * exp(-x); }
auto eqn_4(double x) -> double {
  return 1000.0 * exp(-1 / x) * exp(-1 / (1 - x));
}

int main() {
  std::println("");
  std::println("Equation 1 (simpson): {}", simpson(eqn_1, 0.0, 1.0));
  std::println("Equation 1 (trapezoidal): {}", trapezoidal(eqn_1, 0.0, 1.0));
  std::println("Equation 1 (midpoint): {}", midpoint(eqn_1, 0.0, 1.0));

  std::println("");
  std::println("Equation 2 (simpson): {}", simpson(eqn_2, 0.0, 1.0));
  std::println("Equation 2 (trapezoidal): {}", trapezoidal(eqn_2, 0.0, 1.0));
  std::println("Equation 2 (midpoint): {}", midpoint(eqn_2, 0.0, 1.0));

  std::println("");
  std::println("Equation 3 (simpson): {}", simpson(eqn_3, 0.0, 1.0));
  std::println("Equation 3 (trapezoidal): {}", trapezoidal(eqn_3, 0.0, 1.0));
  std::println("Equation 3 (midpoint): {}", midpoint(eqn_3, 0.0, 1.0));

  std::println("");
  std::println("Equation 4 (simpson): {}", simpson(eqn_4, 0.0, 1.0));
  std::println("Equation 4 (trapezoidal): {}", trapezoidal(eqn_4, 0.0, 1.0));
  std::println("Equation 4 (midpoint): {}", midpoint(eqn_4, 0.0, 1.0));
  return 0;
}
