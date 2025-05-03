#include "nr3.h"
#include "tridag.h"
#include "utilities.h"
#include <print>
#include <vector>

VecDoub target_values(std::vector<VecDoub> y_estimates, const double target_x,
                      const double x_low, const double x_high) {
  const auto size = y_estimates.size();
  VecDoub estimates(size);
  for (size_t i = 0; i < size; i++) {
    const size_t N = y_estimates[i].size() - 1;
    const double h = (x_high - x_low) / (double)N;
    const size_t index = (target_x - x_low) / h;
    estimates[i] = y_estimates[i][index];
  }
  return estimates;
}

// y(target_x) = ?
VecDoub finite_difference_method(const int starting_N, const double target_x,
                                 const double a, const double b,
                                 const double alpha, const double beta,
                                 double F(double y_prime, double y, double x),
                                 double F_y(double y_prime, double y, double x),
                                 double F_y_prime(double y_prime, double y,
                                                  double x)) {
  std::vector<VecDoub> y_estimates;
  int N = starting_N;
  // Initial guess y_i using linear interpolation
  VecDoub y(N + 1);
  for (int i = 0; i < N + 1; i++) {

    y[i] = alpha + (beta - alpha) / (double)N * (double)i;
  }
  y_estimates.push_back(y);

  bool should_stop = false;
  while (!should_stop) {
    const double h = (b - a) / (double)N;
    // Define J (using tridag)
    VecDoub J_a(N + 1, 0.0), J_b(N + 1, 0.0), J_c(N + 1, 0.0);

    {
      // First element
      double x_1 = a + h;
      J_a[0] = 0.0; // Not used
      J_b[0] = 2.0 + pow(h, 2) * F_y((y[2] - alpha) / (2.0 * h), y[1], x_1);
      J_c[0] =
          -1.0 + (h / 2.0) * F_y_prime((y[2] - alpha) / (2.0 * h), y[1], x_1);

      // Last element
      double x_N_1 = b - h;
      J_a[N] = -1. - (h / 2.0) *
                         F_y_prime((beta - y[N - 1]) / (2.0 * h), y[N], x_N_1);
      J_b[N] =
          2.0 + pow(h, 2) * F_y((beta - y[N - 1]) / (2.0 * h), y[N], x_N_1);
      J_c[N] = 0.0; // Not used

      // Middle elements
      for (int i = 1; i < N; i++) {
        double x_i = a + (i + 1) * h;
        J_a[i] =
            -1. -
            h / 2.0 * F_y_prime((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
        J_b[i] =
            2.0 + pow(h, 2) * F_y((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
        J_c[i] =
            -1.0 +
            h / 2.0 * F_y_prime((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
      }
    }

    VecDoub J_a_correct(N - 1), J_b_correct(N - 1), J_c_correct(N - 1),
        J_u_correct(N - 1), J_r_correct(N - 1);
    for (int i = 1; i < N; i++) {
      J_a_correct[i - 1] = J_a[i];
      J_b_correct[i - 1] = J_b[i];
      J_c_correct[i - 1] = J_c[i];
    }
    // Calculate phi(y)=0
    {
      double x_1 = a + h;
      J_r_correct[0] = -alpha + 2.0 * y[1] - y[2] +
                       pow(h, 2) * F((y[2] - alpha) / (2.0 * h), y[1], x_1);
    }
    {
      double x_N_1 = b - h;
      J_r_correct[N - 2] =
          -y[N - 2] - 2.0 * y[N - 1] - beta +
          pow(h, 2) * F((beta - y[N - 2]) / (2.0 * h), y[N - 1], x_N_1);
    }
    for (int i = 2; i < N - 1; i++) {
      double x_i = a + (i + 1) * h;
      J_r_correct[i - 1] =
          -y[i - 1] - 2.0 * y[i] - y[i + 1] +
          pow(h, 2) * F((y[i + 1] - y[i - 1]) / (2.0 * h), y[i], x_i);
    }

    // Use tridag to calculate delta y
    tridag(J_a_correct, J_b_correct, J_c_correct, J_r_correct, J_u_correct);

    for (int i = 1; i < N; i++) {
      y[i] += J_u_correct[i - 1];
    }

    // Half the step size and fill in missing y values with linear
    // interpolation
    N *= 2;
    VecDoub y_new(N + 1);
    y_new[0] = y[0];
    y_new[N] = y[N / 2];
    for (int i = 1; i < N + 1; i++) {
      if (i % 2 == 0) {
        y_new[i] = y[i / 2];
      } else {
        y_new[i] = 0.5 * (y[(i - 1) / 2.0] + y[(i + 1) / 2.0]);
      }
    }

    y = y_new;
    y_estimates.push_back(y);
    if (y_estimates.size() > 5) {
      const auto target_estimates = target_values(y_estimates, target_x, a, b);
      util::print(target_estimates, "target_estimates");
      should_stop = true;
    }
  }

  // Use richardson

  return y;
}

// y_mm(x)  = 2*x+sin(y_m(x))-cos(y(x)) for 0<x<2
// y(0)=0
// y(2)=1

// y''=F
double F(double y_prime, double y, double x) {
  return 2 * x + sin(y_prime) - cos(y);
};

// d/dy F = F_y
double F_y(double y_prime, double y, double x) {
  (void)y_prime; // Unused
  (void)y;       // Unused
  return sin(x);
};

// d/d(y') F = F_y_prime
double F_y_prime(double y_prime, double y, double x) {
  (void)y; // Unused
  (void)x; // Unused
  return cos(y_prime);
};

int main() {
  {
    double a = 0;
    double b = 2;
    double alpha = 0;
    double beta = 1;
    int N = 10;
    // y(1) = ?
    double target_x = 1;
    auto res = finite_difference_method(N, target_x, a, b, alpha, beta, F, F_y,
                                        F_y_prime);
    // util::print(res);
  }
}
