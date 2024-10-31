#include "FFT1DRecursive.H"

#include <cassert>
#include <cmath>

FFT1DRecursive::FFT1DRecursive(const unsigned int &a_M) : FFT1D(a_M) {
  // do you need to build any data structures to help implement this class?
}

FFT1DRecursive::~FFT1DRecursive() {
  // be sure to clean up anything you created in the constructor.
  // most STL containers will clean themselves up for you.
}

void FFT1DRecursive::forwardFFTCC(vector<complex<double>> &a_fHat,
                                  const vector<complex<double>> &f) const {
  assert(f.size() == m_N);
  a_fHat.resize(m_N, 0.0);

  for (int i = 0; i < m_N; i++) {
    a_fHat[i] = forwardRec(f, m_N, i, -1);
  }
}

void FFT1DRecursive::inverseFFTCC(vector<complex<double>> &a_f,
                                  const vector<complex<double>> &a_fHat) const {
  assert(a_fHat.size() == m_N);
  a_f.resize(m_N, 0.0);

  for (int i = 0; i < m_N; i++) {
    a_f[i] = forwardRec(a_fHat, m_N, i, 1);
  }
}

complex<double> FFT1DRecursive::forwardRec(const vector<complex<double>> &f,
                                           int N, int k, int sign) const {
  if (N == 2) {
    if (k == 0)
      return f[0] + f[1];
    else
      return f[0] - f[1];
  }

  vector<complex<double>> even_f(N / 2, 0.0);
  vector<complex<double>> odd_f(N / 2, 0.0);
  for (int i = 0; i < N; i += 2) {
    even_f[i / 2] = f[i];
  }
  for (int i = 1; i < N; i += 2) {
    odd_f[i / 2] = f[i];
  }

  bool pos = true;
  if (k >= N / 2) {
    k -= N / 2;
    pos = false;
  }

  std::complex<double> omega =
      std::exp(std::complex<double>(0, sign * 2.0 * M_PI * k / N));

  complex<double> ans = forwardRec(even_f, N / 2, k, sign);
  if (pos) {
    ans += omega * forwardRec(odd_f, N / 2, k, sign);
  } else {
    ans -= omega * forwardRec(odd_f, N / 2, k, sign);
  }

  return ans;
}
