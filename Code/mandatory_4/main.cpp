#include "nr3.h"
#include "ode.h"
#include "utilities.h"
#include <print>

VecDoub derivs(const Doub t, VecDoub_I &v) {
  VecDoub_O dvdt(3);
  dvdt[0] = exp(-t) * cos(v[1]) + pow(v[2], 2) - v[0];
  dvdt[1] = cos(pow(v[2], 2)) - v[1];
  dvdt[2] = cos(t) * exp(-pow(v[0], 2)) - v[2];
  return dvdt;
}

int main() {
  {
    std::println("\n{:/^80}", " i ");
    double vals[3] = {1.0, 2.0, 3.0};
    double t = 0.0;
    const VecDoub v(3, vals);
    const VecDoub res = derivs(t, v);
    util::print(res, "v1'(0), v2'(0), v3'(0)");
  }
  {
    std::println("\n{:/^80}", " ii ");
    double vals[3] = {1.0, 2.0, 3.0};
    const VecDoub v(3, vals);
    double x_low = 0.0;
    double x_high = 5.0;
    std::vector<VecDoub> result;
    for (int N = 50; N <= 800; N *= 2) {
      auto v_n = v;
      VecDoub res = trapezoidal(x_low, x_high, N, v_n, derivs);
      util::print(res, std::format("Trapezoidal with N={}, h={}", N,
                                   (x_high - x_low) / (double)N));
    }
  }
}
