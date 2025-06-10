#include "nr3.h"
#include "utils/integrals.h"
#include <print>

double eqn(double x) { return exp(pow(x, 3)) * sqrt(x * (2 - x)); }

int main() {
  {
    // 1)
    std::println("\n{:/^100}", " i ");

    std::println("\nSimpson");
    print_quadrature_table(eqn, 0.0, 2.0, IntegrationType::Simpson, 0,
                           pow(2, 20) + 1);
  }

  {
    // 2)
    std::println("\n{:/^100}", " ii ");
  }

  {
    // 3)
    std::println("\n{:/^100}", " iii ");
  }

  {
    // 4)
    std::println("\n{:/^100}", " iv ");
  }

  {
    // 5)
    std::println("\n{:/^100}", " v ");
    std::println("\nTrapezoidal");
    print_quadrature_table(eqn, 0.0, 2.0, IntegrationType::Trapezoidal, 0,
                           pow(2, 20) + 1);
    std::println("\nMidpoint");
    print_quadrature_table(eqn, 0.0, 2.0, IntegrationType::Midpoint, 0,
                           pow(2, 20) + 1);
  }
}
