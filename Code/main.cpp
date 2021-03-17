#include <iostream>
#include "mmsystem.hpp"

int test = 1;

double xl = 0, xr = 2, yl = 0, yr = 2;

double k1(double x, double y)
{
  switch (test) {
    case 1:
      return 1;
    case 2:
      return x + y;
    case 3:
      return x * x + y * y;
    case 4:
      return x * x + y * y;
    default:
      return 1;
  }
}

double k2(double x, double y)
{
  switch (test) {
    case 1:
      return 1;
    case 2:
      return x + y;
    case 3:
      return x * x + y * y;
    case 4:
      return x * x + y * y;
    default:
      return 1;
  }
}

double f(double x, double y)
{
  switch (test) {
    case 1:
      return 0;
    case 2:
      return -2;
    case 3:
      return -8 * x * x - 8 * y * y;
    case 4:
      return -12 * x * x * y * y - 2 * pow(y, 4) - 2 * pow(x, 4);
    default:
      return 1;
  }
}

double u(double x, double y)
{
  switch (test) {
    case 1:
      return 1;
    case 2:
      return x + y;
    case 3:
      return x * x + y * y;
    case 4:
      return x * x * y * y;
    default:
      return 1;
  }
}

double xi1(double y)
{
  switch (test) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    case 4:
      return 4;
    default:
      return 1;
  }
}

double xi4(double x)
{
  switch (test) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    case 4:
      return 4;
    default:
      return 1;
  }
}

double mu1(double y)
{
  switch (test) {
    case 1:
      return xi1(y) * u(xl, y);
    case 2:
      return xi1(y) * u(xl, y) - k1(xl, y);
    case 3:
      return xi1(y) * u(xl, y) - k1(xl, y) * 2 * xl;
    case 4:
      return xi1(y) * u(xl, y) - k1(xl, y) * 2 * xl * y * y;
    default:
      return 1;
  }
}

double mu2(double y)
{
  switch (test) {
    case 1:
      return u(xr, y);
    case 2:
      return u(xr, y);
    case 3:
      return u(xr, y);
    case 4:
      return u(xr, y);
    default:
      return 1;
  }
}

double mu3(double x)
{
  switch (test) {
    case 1:
      return u(x, yl);
    case 2:
      return u(x, yl);
    case 3:
      return u(x, yl);
    case 4:
      return u(x, yl);
    default:
      return 1;
  }
}

double mu4(double x)
{
  switch (test) {
    case 1:
      return xi4(x) * u(x, yr);
    case 2:
      return xi4(x) * u(x, yr) + k2(x, yr);
    case 3:
      return xi4(x) * u(x, yr) + k2(x, yr) * 2 * yr;
    case 4:
      return xi4(x) * u(x, yr) + k2(x, yr) * 2 * yr * x * x;
    default:
      return 1;
  }
}

int main()
{
  bool bools[2] = { true, false };
  MainSystem mainSystem;
  double lastErr[2], currErr[2];

  for (test = 1; test < 5; test++) {
    int Nr = 4, Nz = 4;
    for (int i = 0; i < 7; i++) {
      printf("%d\t%d\t%d\t", Nr, Nz, Nr * Nz);
      for (int b = 0; b < 2; b++) {
        mainSystem.init(xl, xr, yl, yr, Nr, Nz, k1, k2, f, mu1, mu2, mu3, mu4, xi1, xi4);
        mainSystem.printInfo = false;
        mainSystem.build();
        mainSystem.build_L();
        printf("%d\t", mainSystem.solve(10e-7, bools[b]));
        mainSystem.finish();
        currErr[b] = mainSystem.compare(u);
        printf("%10.4e\t", currErr[b]);
        if (test > 2) {
          if (i > 0) printf("%5.4f\t", lastErr[b] / currErr[b]);
          else printf("-\t");
          lastErr[b] = currErr[b];
        }
      }
      printf("\n");
      Nr *= 2;
      Nz *= 2;
    }
    printf("\n\n");
  }
  return 0;
}

