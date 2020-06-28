// zheev.cc
// ʣ�ǥ���ߡ��ȹ���θ�ͭ�����������

#include <stdio.h>
#include <complex>
typedef std::complex<double> Complex;

extern "C" {
  void zheev_ ( const char& JOBZ, const char& UPLO,
		const int& N, Complex** A, const int& LDA,
		double* W, Complex* WORK, const int& LWORK,
		double* RWORK,
		int& INFO, int JOBZlen, int UPLOlen );
};

// ʣ�ǥ���ߡ��ȹ��� H �θ�ͭ�������� H x = E x ��򤯴ʰ״ؿ�
// Input: H[N][N],  Output: U[N][N], E[N]
// ��������H[i][j<=i] �β����Ѥ������Ȥ���ʤ�
// U^\dagger H U ���гѹ���Ȥʤ�
//
template <int N> int zheev( Complex H[N][N], Complex U[N][N], double E[N] )
{
  int i, j, info;
  static Complex cwork[4*N];
  static double  rwork[4*N];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      U[i][j] = H[j][i];
    }
  }

  zheev_( 'V', 'L', N, (Complex**)U, N, E, cwork, 4*N, rwork, info, 1, 1 );

  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      Complex temp;
      temp = U[i][j];
      U[i][j] = U[j][i];
      U[j][i] = temp;
    }
  }

  return info;
}

int main(void)
{
  const int N = 4;
  int i, j, info;

  Complex H[N][N], U[N][N];
  double E[N];

  for( i=0; i<N; i++ ){
    for( j=0; j<=i; j++ ){
      if( i == j+1 ){
	H[i][j] = -1.0;
      }else{
	H[i][j] = 0.0;
      }
    }
  }

  info = zheev( H, U, E );

  printf("# info=%d.\n", info );

  printf("# Eigen values.\n");
  for( j=0; j<N; j++ ){
    printf("%+f ", E[j] );
  }
  printf("\n");

  printf("# Eigen vectors.\n");
  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      printf("%+f%+f ", real(U[i][j]), imag(U[i][j]) );
    }
    printf("\n");
  }

  return 0;
}
