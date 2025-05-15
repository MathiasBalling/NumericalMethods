#include "nr3.h"
#include "utils/ode.h"

VecDoub derivs(const Doub t, VecDoub_I &y) {
  (void)t; // Unused
  VecDoub_O dydt(3);
  dydt[0] = -0.05 * y[0] + 10000.0 * y[1] * y[2];
  dydt[1] = 0.05 * y[0] - 10000.0 * y[1] * y[2] - 30000000.0 * pow(y[1], 2);
  dydt[2] = 30000000.0 * pow(y[1], 2);
  return dydt;
}

int main(int argc, char *argv[]) {

  double vals[3] = {1.0, 0.0, 0.0};
  VecDoub y(3, vals);

  auto res = trapezoidal(0, 5, 1000, y, derivs);

  // Should be around {0.8713, 0.000022, 0.1286}
  util::print(res, "y1(5), y2(5), y3(5)");
}
