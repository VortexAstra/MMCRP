#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>

using namespace std;

class Vect
{
private:
  double * x;
  int n;
public:
  Vect() {}

  Vect(int N)
  {
    n = N;
    x = new double[n];
  }

  Vect(double * X, int N)
  {
    x = X;
    n = N;
  }

  double get(int i)
  {
    return x[i];
  }

  double & operator[](int i)
  {
    return *(x + i);
  }

  void print(FILE * F, char * f)
  {
    for (int i = 0;; i++) {
      fprintf(F, f, x[i]);
      if (i >= n - 1) break;
      fprintf(F, " ");
    }
  }

  double operator*(Vect & q)
  {
    double s = 0;
    for (int i = 0; i < n; i++) s += x[i] * q[i];
    return s;
  }

  Vect clone()
  {
    Vect v(n);
    for (int i = 0; i < n; i++) v[i] = x[i];
    return v;
  }

  Vect & operator=(const Vect & q)
  {
    if (this == &q) {
      return *this;
    }
    n = q.n;
    x = q.x;
    return *this;
  }
};

class Matr
{
private:
  double * m;
public:
  int x, y;

  Matr() {}

  Matr(int x, int y)
  {
    this->x = x;
    this->y = y;
    m = new double[x * y];
  }

  double * operator[](int i) { return (m + y * i); }
};

class Matr1
{
private:
  int nx, ny, n;
public:
  double * C, * B, * A;
  int t;
  Matr1() {}

  Matr1(int x, int y, int type = 0)
  {
    t = type;
    nx = x;
    ny = y;
    n = nx * ny;
    C = new double[n - nx];
    B = new double[n - 1];
    A = new double[n];
  }

  Matr1(double * a, double * b, double * c, int x, int y, int type = 0)
  {
    t = type;
    nx = x;
    ny = y;
    n = nx * ny;
    C = a;
    B = b;
    A = c;
  }

  double get(int i, int j)
  {
    if (abs(i - j) == 0) return A[i];
    else if ((i - j == 1 && t >= 0) || (i - j == -1 && t <= 0)) return B[min(i, j)];
    else if ((i - j == nx && t >= 0) || (i - j == -nx && t <= 0)) return C[min(i, j)];
    else return 0;
  }

  Vect multiplyByVector(Vect x)
  {
    Vect v(n);
    for (int i = 0; i < n; i++) {
      v[i] = A[i] * x[i];
      if (i > 0 && t >= 0) v[i] += B[i - 1] * x[i - 1];
      if (i < n - 1 && t <= 0) v[i] += B[i] * x[i + 1];
      if (i >= nx && t >= 0) v[i] += C[i - nx] * x[i - nx];
      if (i < n - nx && t <= 0) v[i] += C[i] * x[i + nx];
    }
    return v;
  }

  Vect operator*(Vect x)
  {
    return multiplyByVector(x);
  }

  void print(FILE * F, char * f)
  {
    for (int i = 0;; i++) {
      for (int j = 0;; j++) {
        fprintf(F, f, get(i, j));
        if (j >= n - 1) break;
        fprintf(F, " ");
      }
      if (i >= n - 1) break;
      fprintf(F, "\n");
    }
  }
};

class MainSystem
{
private:
public:
  Matr1 A; // Матрица системы
  Matr1 L;
  Vect v; // Решение системы
  Vect g; // Правая часть системы
  Matr u; // Искомые решения
  int Nx, Ny; // Размерность задачи
  int nx, ny, n; // Размерность системы
  double x1, x2, y1, y2, hx, hy; // Границы и шаг
  double (* _k1)(double, double);
  double (* _k2)(double, double);
  double (* _f)(double, double);
  double (* _mu1)(double);
  double (* _mu2)(double);
  double (* _mu3)(double);
  double (* _mu4)(double);
  double (* _xi1)(double);
  double (* _xi4)(double);

  double x(double i) { return x1 + i * hx; }

  double y(double i) { return y1 + i * hy; }

  double f(int i, int j) { return (*_f)(x(i), y(j)); }

  double k1(double i, int j) { return (*_k1)(x(i), y(j)); }

  double k2(int i, double j) { return (*_k2)(x(i), y(j)); }

  double mu1(int j) { return (*_mu1)(y(j)); }

  double mu2(int j) { return (*_mu2)(y(j)); }

  double mu3(int i) { return (*_mu3)(x(i)); }

  double mu4(int i) { return (*_mu4)(x(i)); }

  double xi1(int j) { return (*_xi1)(y(j)); }

  double xi4(int i) { return (*_xi4)(x(i)); }


  bool printInfo = false;

  void init(double xl, double xr, double yl, double yr, int Nx, int Ny,
            double (* k1)(double, double), double (* k2)(double, double),
            double (* f)(double, double), double (* mu1)(double),
            double (* mu2)(double), double (* mu3)(double),
            double (* mu4)(double), double (* xi1)(double), double (* xi4)(double))
  {
    this->Nx = Nx;
    this->Ny = Ny;
    x1 = xl;
    x2 = xr;
    y1 = yl;
    y2 = yr;
    hx = (x2 - x1) / Nx;
    hy = (y2 - y1) / Ny;
    _k1 = k1;
    _k2 = k2;
    _f = f;
    _mu1 = mu1;
    _mu2 = mu2;
    _mu3 = mu3;
    _mu4 = mu4;
    _xi1 = xi1;
    _xi4 = xi4;
    nx = Nx; // -1
    ny = Ny; // -1
    n = nx * ny;
  }

  void build()
  {
    int i, j;
    u = Matr(Nx + 1, Ny + 1);
    g = Vect(n);
    A = Matr1(nx, ny);
    for (j = 1; j <= Ny; j++)
      for (i = 0; i < Nx; i++) {
        int z = nx * (j - 1) + i;
        //i=i; // Одно известное решение при i=Nx
        //j=j-1; // Одно известное решение при j=0
        if (j > 1 && i > 0 && j < ny && i < nx - 1) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy * k1(i - 0.5, j) / hx + hx * k2(i, j + 0.5) / hy + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy * k1(i - 0.5, j) / hx;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy;
          g[z] = hx * hy * f(i, j);
        } else if (j == 1 && i > 0 && i < nx - 1) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy * k1(i - 0.5, j) / hx + hx * k2(i, j + 0.5) / hy + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy * k1(i - 0.5, j) / hx;
          g[z] = hx * hy * f(i, j) + hx * k2(i, j - 0.5) / hy * mu3(i);
        } else if (j == ny && i > 0 && i < nx - 1) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx / 2 + hy * k1(i - 0.5, j) / hx / 2 + hx * xi4(i) + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy * k1(i - 0.5, j) / hx / 2;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy;
          g[z] = hx * hy * f(i, j) / 2 + hx * mu4(i);
        } else if (i == 0 && j > 1 && j < ny) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy *xi1(j) + hx * k2(i, j + 0.5) / hy / 2 + hx * k2(i, j - 0.5) / hy / 2;
          A.B[z - 1] = 0;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy / 2;
          g[z] = hx * hy * f(i, j) / 2 + hy * mu1(j);
        } else if (i == nx - 1 && j > 1 && j < ny) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy * k1(i - 0.5, j) / hx + hx * k2(i, j + 0.5) / hy + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy * k1(i - 0.5, j) / hx;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy;
          g[z] = hx * hy * f(i, j) + hy * k1(i + 0.5, j) / hx * mu2(j);
        } else if (j == 1 && i == 0) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy *xi1(j) + hx * k2(i, j + 0.5) / hy / 2 + hx * k2(i, j - 0.5) / hy / 2;
          g[z] = hx * hy * f(i, j)  / 2 + hy * mu1(j) + hx * k2(i, j - 0.5) / hy * mu3(i) / 2;
        } else if (j == ny && i == 0) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx / 2 + hy * xi1(j) + hx * xi4(i) + hx * k2(i, j - 0.5) / hy / 2;
          A.B[z - 1] = 0;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy / 2;
          g[z] = hx * hy * f(i, j) / 4 + hy * mu1(j) + hx * mu4(i);
        } else if (j == 1 && i == nx - 1) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx + hy * k1(i - 0.5, j) / hx + hx * k2(i, j + 0.5) / hy + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy * k1(i - 0.5, j) / hx;
          g[z] = hx * hy * f(i, j) + hy * k1(i + 0.5, j) / hx * mu2(j) + hx * k2(i, j - 0.5) / hy * mu3(i);
        } else if (j == ny && i == nx - 1) {
          A.A[z] =
          hy * k1(i + 0.5, j) / hx / 2 + hy * k1(i - 0.5, j) / hx / 2 + hx * xi4(i) + hx * k2(i, j - 0.5) / hy;
          A.B[z - 1] = -hy / 2 * k1(i - 0.5, j) / hx;
          A.C[z - nx] = -hx * k2(i, j - 0.5) / hy;
          g[z] = hx / 2 * hy * f(i, j) + hx * mu4(i) + hy * k1(i + 0.5, j) / hx * mu2(j) / 2;
        }
      }
  }

  void build_L()
  {
    int i;
    L = Matr1(nx, ny, 1);
    L.A[0] = sqrt(A.A[0]);
    for (i = 1; i < nx; i++) {
      L.B[i - 1] = A.B[i - 1] / L.A[i - 1];
      L.A[i] = sqrt(A.A[i] - L.B[i - 1] * L.B[i - 1]);
    }
    for (i = nx; i < n; i++) {
      L.C[i - nx] = A.C[i - nx] / L.A[i - nx];
      L.B[i - 1] = A.B[i - 1] / L.A[i - 1];
      L.A[i] = sqrt(A.A[i] - L.B[i - 1] * L.B[i - 1] - L.C[i - nx] * L.C[i - nx]);
    }
  }

  void get_w(Vect r, Vect & w)
  {
    int i;
    Vect y(n);
    y[0] = r[0] / L.A[0];
    for (i = 1; i < nx; i++) y[i] = (r[i] - L.B[i - 1] * y[i - 1]) / L.A[i];
    for (i = nx; i < n; i++) y[i] = (r[i] - L.B[i - 1] * y[i - 1] - L.C[i - nx] * y[i - nx]) / L.A[i];
    w[n - 1] = y[n - 1] / L.A[i - 1];
    for (i = n - 2; i >= n - nx; i--) w[i] = (y[i] - L.B[i] * w[i + 1]) / L.A[i];
    for (i = n - nx - 1; i >= 0; i--) w[i] = (y[i] - L.B[i] * w[i + 1] - L.C[i] * w[i + nx]) / L.A[i];
  }

  int solve(double e, bool yav)
  {
    Vect ax, as;
    Vect s(n), x(n), r(n), w(n);
    double a, b, nr, wr;
    build_L();
    for (int i = 0; i < n; i++) x[i] = 0;
    ax = A * x;
    for (int i = 0; i < n; i++) r[i] = g[i] - ax[i];
    if (!yav) get_w(r, w); else w = r;
    for (int i = 0; i < n; i++) s[i] = w[i];
    int l = 0;
    if (r * r != 0)
      do {
        as = A * s;
        wr = w * r;
        a = wr / (as * s);
        for (int i = 0; i < n; i++) {
          x[i] = x[i] + s[i] * a;
          r[i] = r[i] - as[i] * a;
        }
        nr = sqrt(r * r) / sqrt(g * g);
        if (!yav) get_w(r, w); else w = r;
        b = (w * r) / wr;
        for (int i = 0; i < n; i++) s[i] = w[i] + b * s[i];
        l++;
        if (printInfo) printf("%d: %g\n", l, nr);
      } while (nr >= e);
    v = x;
    return l;
  }

  double check()
  {
    double s = 0;
    Vect ax;
    ax = A * v;
    for (int i = 0; i < n; i++) if (s < abs(ax[i] - g[i])) s = abs(ax[i] - g[i]);
    return s;
  }

  void finish()
  {
    for (int i0 = 0; i0 < Nx; i0++) {
      for (int j0 = 1; j0 <= Ny; j0++) {
        u[i0][j0] = v[nx * (j0 - 1) + i0];
        if (printInfo) printf("%9.6f ", u[i0][j0]);
      }
      if (printInfo) printf("\n");
    }
  }

  double check_matrix(double (* u0)(double, double))
  {
    double s = 0;
    Vect au;
    Vect u(n);
    for (int i0 = 0; i0 < Nx; i0++)
      for (int j0 = 1; j0 <= Ny; j0++)
        u[nx * (j0 - 1) + i0] = u0(x(i0), y(j0));
    au = A * u;
    for (int i0 = 0; i0 < Nx; i0++) {
      for (int j0 = 1; j0 <= Ny; j0++) {
        int i = nx * (j0 - 1) + i0;
        if (abs(g[i] - au[i]) > s) s = abs(g[i] - au[i]);
        if (printInfo) printf("%9.6f ", abs(g[i] - au[i]));
      }
      if (printInfo) printf("\n");
    }
    return s;
  }

  double compare(double (* u0)(double, double))
  {
    double s = 0;
    for (int i = 0; i < Nx; i++) {
      for (int j = 1; j <= Ny; j++) {
        s = fmax(s, abs(u0(x(i), y(j)) - u[i][j]));
        if (printInfo) printf("%6.3f ", abs(u0(x(i), y(j)) - u[i][j]));
      }
      if (printInfo) printf("\n");
    }
    return s;
  }

  double get(double xx, double yy)
  {
    double di = (xx - x1) / hx;
    double dj = (yy - y1) / hy;
    int i = (int) di;
    int j = (int) dj;
    return u[i][j];
  }
};


