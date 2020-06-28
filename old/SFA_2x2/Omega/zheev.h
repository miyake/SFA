// zheev.h

#include <complex>
 typedef std::complex<double> Complex;

extern "C" {
  void zheev_ ( const char& JOBZ, const char& UPLO,
		const int& N, Complex** A, const int& LDA,
		double* W, Complex* WORK, const int& LWORK,
		double* RWORK,
		int& INFO);
};




