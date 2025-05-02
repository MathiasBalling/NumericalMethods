#include "de_rule_table.h"
#include "quadrature_table.h"
#include <print>

double eqn(double x) { return (cos(pow(x, 3)) * exp(-x)) / sqrt(x); }

double eqn_derule(double x, double delta) {
  // If x is small use delta instead to avoid division by zero
  // TODO: How to deal with delta
  if (abs(x) < 1e-6) {
    return (cos(pow(x, 3)) * exp(-x)) / sqrt(delta);
  } else {
    return (cos(pow(x, 3)) * exp(-x)) / sqrt(x);
  }
}
double limit_low = 0.0;
double limit_high = 4.0;

int main() {
  {
    // Correct answer: ~1.43004
    std::println("\n{:/^120}", " Extended Midpoint: ");
    print_quadrature_table(eqn, limit_low, limit_high,
                           IntegrationType::Midpoint, 1e-3);
  }

  {
    // Correct answer: ~1.43004
    std::println("\n{:/^120}", " DErule: ");
    auto derule = DEruleTable(eqn_derule, limit_low, limit_high);
    std::println("");
    double error = std::numeric_limits<double>::max();
    double diff = std::numeric_limits<double>::max();
    double prev = 0;
    while (abs(diff) > 1e-03) {
      double res = derule.next();
      error = derule.get_error();
      diff = res - prev;
      prev = res;
    }
    derule.print_table();
  }
}
