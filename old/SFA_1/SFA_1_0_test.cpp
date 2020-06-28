// SFA, 1サイト, T>0

#include <iostream.h>
#include <math.h>
#include <process.h>
#include "dsyev.h"
#include <complex>
 typedef std::complex<double> Complex;


main()
{

// 準備
int i,j,k,l,m,n,N,n_omega,z,Npur[4],Npdr[4],Npu[4],Npd[4];
double cu[4][4],cd[4][4],w[4][4],v[4][4];
double e,t,U,mu,ee;
double E[4];
double E0,Omega_r,Rr[4],Rg[4],Omega_t[4],Omega, omega,a; double au[4][4],ad[4][4],bu[4][4],bd[4][4];
double GRu[4],GRd[4],SEu[4],SEd[4]; double GIu[4][64][64],GId[4][64][64],gIu[64],gId[64];
double Fu[64][64],Fd[64][64],work[191]; char jobz='V'; char uplo='U'; int info;
// 生成消滅演算子、基底ベクトル、基底状態
for(i=0;i<=3;++i){
	for(j=0;j<=3;++j){
		cu[i][j]=cd[i][j]=w[i][j]=v[i][j]=0.0;
	};
};
cu[0][1]=cu[2][3]=cd[0][2]=1.0; cd[1][3]=-1.0;
for(i=0;i<=3;++i){
	w[i][i]=1.0;
};
for(i=0;i<=3;++i){
	v[0][i]=w[1][i]; v[1][i]=w[2][i];
	v[2][i]=w[1][i]/(sqrt(2))+w[2][i]/(sqrt(2));
	v[3][i]=w[1][i]/(sqrt(2))-w[2][i]/(sqrt(2));
};
// 代理系のグリーン関数の重み
for(z=0;z<=3;++z){
	for(n=0;n<=3;++n){
		for(i=0;i<=3;++i){
			for(j=0;j<=3;++j){
				au[z][n]+=v[z][i]*cu[i][j]*w[n][j]; bu[z][n]+=w[n][j]*cu[j][i]*v[z][i];
				ad[z][n]+=v[z][i]*cd[i][j]*w[n][j]; bd[z][n]+=w[n][j]*cd[j][i]*v[z][i];
			};
		};
	};
};
// 元の系のパラメーター
e=0.0; t=1.0;
cout << " (0<mu<U) の範囲で\n";
cout << "U="; cin >> U;
cout << "mu="; cin >> mu;
// 代理系のパラメーター
ee=0.0;
// 代理系のハミルトニアン、代理系の熱力学ポテンシャル
E[0]=0.0; E[1]=ee-mu; E[2]=ee-mu; E[3]=2*ee+U-2*mu;
E0=Omega_r=ee-mu;
// z スタート
for(z=0;z<=3;++z){
cout << "-------(z="<<z<<")--------\n";
Npur[z]=Npdr[z]=Npu[z]=Npd[z]=0;
// omega スタート
for(n_omega=0;n_omega<=5000;++n_omega){
  omega=-5.0+0.001*n_omega;
  // 代理系のグリーン関数
  GRu[z]=GRd[z]=0.0;
  for(n=0;n<=3;++n){
  	if(fabs(au[z][n])>0.0000001){
  		GRu[z]+=au[z][n]*au[z][n]/(omega+E0-E[n]); Rr[z]+=-E0+E[n]; Npur[z]+=1;
  	};
  	if(fabs(bu[z][n])>0.0000001){
  		GRu[z]+=bu[z][n]*bu[z][n]/(omega+E[n]-E0); Rr[z]+=-E[n]+E0; Npur[z]+=1;
  	};
  	if(fabs(ad[z][n])>0.0000001){
  		GRd[z]+=ad[z][n]*ad[z][n]/(omega+E0-E[n]); Rr[z]+=-E0+E[n]; Npdr[z]+=1;
  	};
  	if(fabs(bd[z][n])>0.0000001){
  		GRd[z]+=bd[z][n]*bd[z][n]/(omega+E[n]-E0); Rr[z]+=-E[n]+E0; Npdr[z]+=1;
  	};
  };
  // 代理系の自己エネルギー
  SEu[z]=omega+mu-ee-GRu[z]; SEd[z]=omega+mu-ee-GRd[z];
  // 元の系のグリーン関数
  for(i=0;i<=63;++i){
	for(j=0;j<=63;++j){
		GIu[z][i][j]=GId[z][i][j]=0.0;
	};
  };
  for(i=0;i<=62;++i){
  	 	GIu[z][i][i]=omega+mu-e-SEu[z]; GId[z][i][i]=omega+mu-e-SEd[z];
  	 	GIu[z][i][i+1]=GIu[z][i+1][i]=GId[z][i][i+1]=GIu[z][i+1][i]=t;
  };
  GIu[z][63][63]=omega+mu-e-SEu[z]; GId[z][63][63]=omega+mu-e-SEd[z];
  GIu[z][0][63]=GIu[z][63][0]=GId[z][0][63]=GId[z][63][0]=t;
  for(i=0;i<=63;++i){
	for(j=0;j<=63;++j){
		Fu[i][j]=GIu[z][i][j]; Fd[i][j]=GId[z][i][j];
	};
  };
  dsyev_( jobz, uplo, 64, (double**)Fu, 64, gIu, work, 191, info);
  if(info!=0) cout << "info=" << info << " (gIu:z="<<z<<", omega="<<omega<<")\n";
  dsyev_( jobz, uplo, 64, (double**)Fd, 64, gId, work, 191, info);
  if(info!=0) cout << "info=" << info << " (gId:z="<<z<<", omega="<<omega<<")\n";
  // Rg
  for(i=0;i<=63;++i){
  	if(fabs(gIu[i])<0.0005){
		Rg[z]+=omega/64.0; Npu[z]+=1;
	};
  	if(fabs(gId[i])<0.0005){
		Rg[z]+=omega/64.0; Npd[z]+=1;
	};
  };
}; // omega のループ、終わり
// Omega_t
Omega_t[z]=Omega_r-Rr[z]+Rg[z];
cout <<"Omega_t[z]="<<Omega_t[z]<<" ( Omega_r="<<Omega_r<<", Rr[z]="<<Rr[z]<<", Rg[z]="<<Rg[z]<<" )\n";
cout <<"Npur="<<Npur[z]<<", Npdr="<<Npdr[z]<<", Npu="<<Npu[z]<<", Npd="<<Npd[z]<<"\n";
}; // z のループ、終わり

Complex b;
b=(1.0,1.0);
cout << b <<"\n";




};
