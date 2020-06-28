// dsyev.cc
// 実対称行列の固有値方程式を解く

#include <stdio.h>

extern "C" {
  void dsyev_ ( const char& JOBZ, const char& UPLO,
		const int& N, double** A, const int& LDA,
		double* W, double* WORK, const int& LWORK,
		int& INFO, int JOBZlen, int UPLOlen );
};

// 実対称行列 H の固有値方程式 H x = E x を解く簡易関数
// Input: H[N][N],  Output: U[N][N], E[N]
// ただし、H[i][j<=i] の下三角しか参照されない
// U^\dagger H U が対角行列となる
//
template <int N> int dsyev( double H[N][N], double U[N][N], double E[N] )
{
  int i, j, info;
  static double work[4*N];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      U[i][j] = H[j][i];
    }
  }

  dsyev_( 'V', 'L', N, (double**)U, N, E, work, 4*N, info, 1, 1 );

  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      double temp;
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

  double H[N][N], U[N][N], E[N];

  for( i=0; i<N; i++ ){
    for( j=0; j<=i; j++ ){
      if( i == j+1 ){
	H[i][j] = -1.0;
      }else{
	H[i][j] = 0.0;
      }
    }
  }

  info = dsyev( H, U, E );

  printf("# info=%d.\n", info );

  printf("# Eigen values.\n");
  for( j=0; j<N; j++ ){
    printf("%+f ", E[j] );
  }
  printf("\n");

  printf("# Eigen vectors.\n");
  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      printf("%+f ", U[i][j] );
    }
    printf("\n");
  }

  return 0;
}
