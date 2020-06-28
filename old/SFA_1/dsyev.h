// dsyev.h

extern "C" {
  void dsyev_ ( const char& JOBZ, const char& UPLO,
		const int& N, double** A, const int& LDA,
		double* W, double* WORK, const int& LWORK,
		int& INFO);
};




