
#include "FFT1DW.H"
#include "fftw3.h"
#include <cstdlib>

FFT1DW::FFT1DW(unsigned int a_N) : FFT1D(a_N) {
  // build any data holders and fftw plans you will need here
}

FFT1DW::~FFT1DW() {
  // clean up any objects you created.  don't forget fftw
  // has it's own clean up routines.
}

void FFT1DW::forwardFFTCC(vector<complex<double>> &a_fHat,
                          const vector<complex<double>> &f) const {
  a_fHat.resize(m_N, 0.0);

  fftw_plan forward;
  fftw_complex *in, *out;
  vector<complex<double>> vec_in = f;

  in = reinterpret_cast<fftw_complex *>(&vec_in[0]);
  out = reinterpret_cast<fftw_complex *>(&a_fHat[0]);

  forward = fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(forward);
  fftw_destroy_plan(forward);
}

void FFT1DW::inverseFFTCC(vector<complex<double>> &a_f,
                          const vector<complex<double>> &a_fHat) const {
  a_f.resize(m_N, 0.0);

  fftw_plan inverse;
  fftw_complex *in, *out;
  vector<complex<double>> vec_in = a_fHat;

  in = reinterpret_cast<fftw_complex *>(&vec_in[0]);
  out = reinterpret_cast<fftw_complex *>(&a_f[0]);

  inverse = fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(inverse);
  fftw_destroy_plan(inverse);
}
