#include "quadrature_table.h"
#include <print>

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
