#ifndef DE_RULE_TABLE_H
#define DE_RULE_TABLE_H
#include "nr3.h"
#include "quadrature.h"
#include <numbers>
template <class T> struct DEruleTable : Quadrature {

  DEruleTable(T &funcc, const Doub aa, const Doub bb, const Doub hmaxx = 3.7)
      : func(funcc), a(aa), b(bb), hmax(hmaxx) {
    n = 0;
    m_error = std::numeric_limits<double>::max();
  }

  Doub next() {
    Doub del, fact, q, sum, t, twoh;
    Int it, j;
    n++;
    if (n == 1) {
      fact = 0.25;
      s = hmax * 2.0 * (b - a) * fact * func(0.5 * (b + a), 0.5 * (b - a));
      f_calls.push_back(1);   // For table
      s_vals.push_back(s);    // For table
      h_vals.push_back(hmax); // For table
      calculate_error();
      return s;
    } else {
      for (it = 1, j = 1; j < n - 1; j++)
        it <<= 1;
      twoh = hmax / it;
      h_vals.push_back(twoh); // For table
      t = 0.5 * twoh;
      int f_sum = f_calls.at(n - 2); // For table
      for (sum = 0.0, j = 0; j < it; j++) {
        q = exp(-2.0 * sinh(t));
        del = (b - a) * q / (1.0 + q);
        fact = q / SQR(1.0 + q) * cosh(t);
        sum += fact * (func(a + del, del) + func(b - del, del));
        f_sum += 2; // For tables
        t += twoh;
      }
      f_calls.push_back(f_sum); // For table
      s = 0.5 * s + (b - a) * twoh * sum;
      s_vals.push_back(s); // For table
      calculate_error();
      return s;
    }
  }

  void print_table() {
    // Print the table
    // Table header
    std::println("|{:^10}|{:^21}|{:^21}|{:^21}|{:^21}|", "i", "A", "A-diff",
                 "e_d", "f-calls");
    // Table header separator
    std::println("|{:-^10}|{:-^21}|{:-^21}|{:-^21}|{:-^21}|", "", "", "", "",
                 "", "");

    // Table body
    for (size_t i = 0; i < s_vals.size(); i++) {
      if (i == 0) {
        std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21}|", i + 1,
                     s_vals.at(i), "*",
                     exp(-pow(std::numbers::pi, 2) / h_vals.at(i)),
                     f_calls.at(i));
      } else {
        std::println("|{:^10}|{:^21.12}|{:^21.12}|{:^21.12}|{:^21}|", i + 1,
                     s_vals.at(i), s_vals.at(i) - s_vals.at(i - 1),
                     exp(-pow(std::numbers::pi, 2) / h_vals.at(i)),
                     f_calls.at(i));
      }
    }
  }

  Doub a, b, hmax, s;
  T &func;
  // Current error estimate
  double get_error() { return m_error; }

private:
  Doub m_error;
  std::vector<int> f_calls;   // Function calls
  std::vector<double> h_vals; // Step sizes
  std::vector<double> s_vals; // Integral values
  //
  void calculate_error() {
    m_error = exp(-pow(std::numbers::pi, 2) / h_vals.back());
  }
};
#endif // !DE_RULE_TABLE_H
