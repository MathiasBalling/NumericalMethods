#ifndef ODE_H
#define ODE_H

#include "multi_roots_table.h"
#include "nr3.h"
#include "utilities.h"

// Euler method
VecDoub
first_order_runge_kuttea_method(double low, double high, int steps, VecDoub_I y,
                                VecDoub derivs(const Doub x, VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;

  VecDoub y_n = y;
  for (double x = low; x < high; x += h) {
    auto y_n_next = y_n + h * derivs(x, y_n);
    y_n = y_n_next;
  }
  return y_n;
}

// Midpoint method
VecDoub second_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &y,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto y_n_next = y_n + k2;
    y_n = y_n_next;
  }
  return y_n;
}

VecDoub fourth_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &y,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = low; x_n < high; x_n += h) {
    auto k1 = h * derivs(x_n, y_n);
    auto k2 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k1);
    auto k3 = h * derivs(x_n + 0.5 * h, y_n + 0.5 * k2);
    auto k4 = h * derivs(x_n + h, y_n + k2);
    auto y_n_next = y_n + 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
    y_n = y_n_next;
  }
  return y_n;
}

VecDoub trapezoidal(double x_low, double x_high, int steps, VecDoub_I &y,
                    VecDoub derivs(const Doub x, VecDoub_I &y)) {
  // Trapezoidal method
  const double h = (x_high - x_low) / (double)steps;
  VecDoub y_n = y;
  for (double x_n = x_low; x_n < x_high; x_n += h) {
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
#endif // !ODE_H
