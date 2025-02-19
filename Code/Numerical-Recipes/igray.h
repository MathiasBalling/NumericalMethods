#ifndef _IGRAY_H_
#define _IGRAY_H_
#include "nr3.h"
struct Gray {

  Uint gray(const Uint n) { return n ^ (n >> 1); }

  Uint invgray(const Uint n) {
    Int ish = 1;
    Uint ans = n, idiv;
    for (;;) {
      ans ^= (idiv = ans >> ish);
      if (idiv <= 1 || ish == 16)
        return ans;
      ish <<= 1;
    }
  }
};
#endif // _IGRAY_H_
