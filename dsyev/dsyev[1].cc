// dsyev.cc
// ���оι���θ�ͭ�����������

#include <stdio.h>

extern "C" {
  void dsyev_ ( const char& JOBZ, const char& UPLO,
		const int& N, double** A, const int& LDA,
		double* W, double* WORK, const int& LWORK,
		int& INFO, int JOBZlen, int UPLOlen );
};

// ���оι��� H �θ�ͭ�������� H x = E x ��򤯴ʰ״ؿ�
// Input: H[N][N],  Output: U[N][N], E[N]
// ��������H[i][j<=i] �β����Ѥ������Ȥ���ʤ�
// U^\dagger H U ���гѹ���Ȥʤ�
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
