#include "quadrature_table.h"
#include <print>

double eqn_1(double x) { return cos(x * x) * exp(-x); }
double eqn_2(double x) { return sqrt(x) * cos(x * x) * exp(-x); }
double eqn_3(double x) { return (1 / sqrt(x)) * cos(x * x) * exp(-x); }
double eqn_4(double x) { return 1000.0 * exp(-1 / x) * exp(-1 / (1 - x)); }

int main() {
  std::println("\nSimpson eqn 1:");
  print_quadrature_table(eqn_1, 0.0, 1.0, IntegrationType::Simpson);

  std::println("\nMidpoint eqn 1:");
  print_quadrature_table(eqn_2, 0.0, 1.0, IntegrationType::Midpoint);

  std::println("\nTrapezoidal eqn 1:");
  print_quadrature_table(eqn_1, 0.0, 1.0, IntegrationType::Trapezoidal);

  std::println("\nSimpson eqn 2:");
  print_quadrature_table(eqn_2, 0.0, 1.0, IntegrationType::Simpson);

  std::println("\nMidpoint eqn 3:");
  print_quadrature_table(eqn_3, 0.0, 1.0, IntegrationType::Midpoint);

  std::println("\nTrapezoidal eqn 4:");
  print_quadrature_table(eqn_4, 0.0, 1.0, IntegrationType::Trapezoidal);
  return 0;
}
