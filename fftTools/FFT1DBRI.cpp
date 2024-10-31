#include "FFT1DBRI.H"
#include "FFT1D.H"
#include "PowerItoI.H"
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <vector>

constexpr double M_PI = 3.14159265358979323846;

FFT1DBRI::FFT1DBRI(const unsigned int &a_M) : FFT1D(a_M) {
  int levelLength = 2;
  int sign;
  int N = Power(2, a_M);
  // if (a_isForward) {sign = -1;} else {sign=1;}
  m_zToTheKF.resize(a_M);
  m_zToTheKI.resize(a_M);
  for (int level = 1; level < a_M; level++) {
    double h = 1. / levelLength / 2;
    m_zToTheKF[level].push_back(complex<double>(1., 0.));
    m_zToTheKI[level].push_back(complex<double>(1., 0.));
    for (int j = 1; j < levelLength; j++) {
      auto zf = exp(complex<double>(0, -2 * M_PI * h * j));
      auto zi = exp(complex<double>(0, 2 * M_PI * h * j));
      m_zToTheKF[level].push_back(zf);
      m_zToTheKI[level].push_back(zi);
    }
    levelLength *= 2;
  }
};
FFT1DBRI::~FFT1DBRI(){};
void FFT1DBRI::forwardFFTCC(vector<complex<double>> &a_fHat,
                            const vector<complex<double>> &a_f) const {
  assert(a_fHat.size() == m_N);
  for (unsigned int j = 0; j < m_N; j++) {
    a_fHat[j] = a_f[j];
  }
  bool isForward = true;

  FFTCTBRI(a_fHat, m_zToTheKF, m_M);
}
void FFT1DBRI::inverseFFTCC(vector<complex<double>> &a_f,
                            const vector<complex<double>> &a_fHat) const {
  assert(a_fHat.size() == m_N);
  for (unsigned int j = 0; j < m_N; j++) {
    a_f[j] = a_fHat[j];
  }
  bool isForward = false;
  FFTCTBRI(a_f, m_zToTheKI, m_M);
}
void FFT1DBRI::FFTCTBRI(vector<complex<double>> &a_f,
                        const vector<vector<complex<double>>> &a_zToTheK,
                        const int &a_M) const {
  int N = Power(2, a_M);

  vector<complex<double>> f(N);

  for (int j = 0; j < N; j++) {
    int jlevel = j;
    int factor = N / 2;
    int indOfj;
    indOfj = 0;
    for (int level = 0; level < a_M; level++) {
      int i = jlevel % 2;
      jlevel = jlevel / 2;
      indOfj += factor * i;
      factor = factor / 2;
    }
    f[indOfj] = a_f[j];
  }
  for (int j = 0; j < N; j += 2) {
    complex<double> fHat0 = f[j] + f[j + 1];
    complex<double> fHat1 = f[j] - f[j + 1];
    f[j] = fHat0;
    f[j + 1] = fHat1;
  }
  int levelLength = 2;
  int sign;
  for (int level = 1; level < a_M; level++) {
    double h = 1. / levelLength / 2;
    for (int j0 = 0; j0 < N; j0 += 2 * levelLength) {
      for (int j = 0; j < levelLength; j++) {
        complex<double> temp = a_zToTheK[level][j] * f[j + j0 + levelLength];
        complex<double> fHat0 = f[j + j0] + temp;
        complex<double> fHat1 = f[j + j0] - temp;
        f[j + j0] = fHat0;
        f[j + j0 + levelLength] = fHat1;
      }
    }
    levelLength *= 2;
  }
  for (int j = 0; j < N; j++) {
    a_f[j] = f[j];
  }
};
