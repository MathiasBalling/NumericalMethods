#include "de_rule_table.h"
#include <print>

double eqn_1(double x, double delta) {
  (void)delta; // unused (No singularity)

  return cos(x * x) * exp(-x);
}

double eqn_2(double x, double delta) {
  (void)delta; // unused (No singularity)
  return sqrt(x) * cos(x * x) * exp(-x);
}

double eqn_3(double x, double delta) {
  if (abs(x) < 1e-10) {
    return (1 / sqrt(delta)) * cos(x * x) * exp(-x);
  } else {
    return (1 / sqrt(x)) * cos(x * x) * exp(-x);
  }
}

double eqn_4(double x, double delta) {
  if (abs(x) < 1e-10) {
    return 1000.0 * exp(-1 / delta) * exp(-1 / (1 - x));
  } else if ((1 - abs(x)) < 1e-10) {
    return 1000.0 * exp(-1 / x) * exp(-1 / (1 - delta));
  } else {
    return 1000.0 * exp(-1 / x) * exp(-1 / (1 - x));
  }
}

int main(int argc, char *argv[]) {
  {
    auto derule = DEruleTable(eqn_1, 0.0, 1.0);
    println("");
    double diff = std::numeric_limits<double>::max();
    double prev = 0;
    while (abs(diff) > 1e-10) {
      double cur = derule.next();
      diff = cur - prev;
      prev = cur;
    }
    derule.print_table();
  }

  {
    auto derule = DEruleTable(eqn_2, 0.0, 1.0);
    println("");
    double diff = std::numeric_limits<double>::max();
    double prev = 0;
    while (abs(diff) > 1e-10) {
      double cur = derule.next();
      diff = cur - prev;
      prev = cur;
    }
    derule.print_table();
  }

  {
    auto derule = DEruleTable(eqn_3, 0.0, 1.0, 4.3);
    println("");
    double diff = std::numeric_limits<double>::max();
    double prev = 0;
    while (abs(diff) > 1e-10) {
      double cur = derule.next();
      diff = cur - prev;
      prev = cur;
    }
    derule.print_table();
  }

  {
    auto derule = DEruleTable(eqn_4, 0.0, 1.0);
    println("");
    double diff = std::numeric_limits<double>::max();
    double prev = 0;
    while (abs(diff) > 1e-10) {
      double cur = derule.next();
      diff = cur - prev;
      prev = cur;
    }
    derule.print_table();
  }
}
