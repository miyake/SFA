// zgesv.cc
// 複素一般行列の線形方程式を解く

#include <stdio.h>
#include <complex>
typedef std::complex<double> Complex;

extern "C" {
  void zgesv_ ( const int& N, const int& NRHS,
		Complex** A, const int& LDA, int* IPIV,
		Complex** B, const int& LDB, int& INFO );
};

// 複素一般行列 A の線形方程式 A x = b を解く簡易関数
// Input: A[N][N], b[N]  Output: x[N]
//
template <int N> int zgesv( Complex A[N][N], Complex x[N], const Complex b[N] )
{
  int i, j, info;
  static int ipiv[N];
  static Complex U[N][N];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      U[i][j] = A[j][i];
    }
    x[i] = b[i];
  }

  zgesv_( N, 1, (Complex**)U, N, ipiv, (Complex**)x, N, info );

  return info;
}

int main(void)
{
  const int N = 4;
  int i, j, info;

  Complex A[N][N];
  Complex x[N], b[N];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      if( i==j ){
	A[i][j] = 2.0;
      }else{
	A[i][j] = 1.0;
      }
    }
    b[i] = 5.0;
  }

  info = zgesv( A, x, b );

  printf("# info=%d.\n", info );

  printf("# Solution.\n");
  for( j=0; j<N; j++ ){
    printf("%+f%+f\n", real(x[j]), imag(x[j]) );
  }

  printf("# Error.\n");
  for( i=0; i<N; i++ ){
    Complex sum=0.0;
    for( j=0; j<N; j++ ){
      sum += A[i][j]*x[j];
    }
    sum -= b[i];

    printf("%+f%+f\n", real(sum), imag(sum) );
  }

  return 0;
}
