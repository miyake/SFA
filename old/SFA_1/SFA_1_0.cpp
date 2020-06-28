// SFA, 1サイト, T>0

#include <iostream.h>
#include <math.h>
#include <process.h>
#include "dsyev.h"

main()
{

int i,j,k,l,m,n,N,n_omega,Npr;
// 生成消滅演算子
double cu[4][4],cd[4][4];
for(i=0;i<=3;++i){
	for(j=0;j<=3;++j){
		cu[i][j]=cd[i][j]=0.0;
	};
};
cu[0][1]=cu[2][3]=cd[0][2]=1.0; cd[1][3]=-1.0;
// 準備
double e,t,U,h, mu, ee,hh,P;
double H[4][4]; double E[4],work[191]; char jobz='V'; char uplo='U'; int info;
double Omega_r,Rr,Rg,Omega_t,Omega, omega,a; double au[4],ad[4],bu[4],bd[4];
double GRu,GRd,SEu,SEd; double GIu[64][64],GId[64][64],gIu[64],gId[64];
// 元の系のパラメーター
e=0.0; t=1.0; h=0.0;
cout << "U="; cin >> U;
cout << "mu="; cin >> mu;
// 代理系のパラメーター
ee=0.0; hh=0.0; P=0.0;
// 初期化
for(i=0;i<=3;++i){
	for(j=0;j<=3;++j){
		H[i][j]=0.0;
	};
	au[i]=ad[i]=bu[i]=bd[i]=0.0;
};
Omega_r=Rr=Rg=0.0; Npr=0;
// 代理系のハミルトニアン
H[0][0]=0.0; H[1][1]=ee-hh-mu; H[2][2]=ee+hh-mu; H[3][3]=2*ee+U-2*mu;
H[0][3]=H[3][0]=P;
// 代理系のハミルトニアンの対角化
dsyev_( jobz, uplo, 4, (double**)H, 4, E, work, 11, info);
if(info!=0) cout << "info=" << info << "\n";
if(E[0]==E[1]){
	cout<<"代理系の基底状態は縮退しています : H[0][i]=("<<
	H[0][0]<<","<<H[0][1]<<","<<H[0][2]<<","<<H[0][3]<<")\n";
};
// 代理系の熱力学ポテンシャル
Omega_r=E[0];
// 生成消滅演算子の基底のとりかえ
for(n=0;n<=3;++n){
	for(i=0;i<=3;++i){
		for(j=0;j<=3;++j){
			au[n]+=H[0][i]*cu[i][j]*H[n][j]; bu[n]+=H[n][j]*cu[j][i]*H[0][i];
			ad[n]+=H[0][i]*cd[i][j]*H[n][j]; bd[n]+=H[n][j]*cd[j][i]*H[0][i];
		};
	};
};
// omega スタート
for(n_omega=0;n_omega<=500;++n_omega){
  omega=-5.0+0.01*n_omega;
  // 代理系のグリーン関数
  GRu=GRd=0.0;
  for(n=0;n<=3;++n){
  	if(fabs(au[n])>0.0000001){
  		GRu+=au[n]*au[n]/(omega+E[0]-E[n]); Rr+=-E[0]+E[n]; Npr+=1;
  	};
  	if(fabs(bu[n])>0.0000001){
  		GRu+=bu[n]*bu[n]/(omega+E[n]-E[0]); Rr+=-E[n]+E[0]; Npr+=1;
  	};
  	if(fabs(ad[n])>0.0000001){
  		GRd+=ad[n]*ad[n]/(omega+E[0]-E[n]); Rr+=-E[0]+E[n]; Npr+=1;
  	};
  	if(fabs(bd[n])>0.0000001){
  		GRd+=bd[n]*bd[n]/(omega+E[n]-E[0]); Rr+=-E[n]+E[0]; Npr+=1;
  	};
  };
  // 代理系の自己エネルギー
  SEu=omega+mu-ee+hh-GRu; SEd=omega+mu-ee+hh-GRd;
  for(i=0;i<=62;++i){
  	 	GIu[i][i]=omega+mu-e+h-SEu; GId[i][i]=omega+mu-e-h-SEd;
  	 	GIu[i][i+1]=GIu[i+1][i]=GId[i][i+1]=GIu[i+1][i]=t;
  };
  GIu[63][63]=omega+mu-e+h-SEu; GId[63][63]=omega+mu-e-h-SEd;
  GIu[0][63]=GIu[63][0]=GId[0][63]=GId[63][0]=t;
  dsyev_( jobz, uplo, 64, (double**)GIu, 64, gIu, work, 191, info);
  if(info!=0) cout << "info=" << info << " (gIu:omega="<<omega<<")\n";
  dsyev_( jobz, uplo, 64, (double**)GId, 64, gId, work, 191, info);
  if(info!=0) cout << "info=" << info << " (gId:omega="<<omega<<")\n";
  // Rg
  for(i=0;i<=63;++i){
  	if(fabs(gIu[i])<0.005) Rg+=omega/64.0;
  	if(fabs(gId[i])<0.005) Rg+=omega/64.0;
  };
}; // omega のループ、終わり
// Omega_t
Omega_t=Omega_r-Rr+Rg;
cout <<"Omega_t="<<Omega_t<<" ( Omega_r="<<Omega_r<<", Rr="<<Rr<<", Rg=" << Rg << " )\n";
  
  
  
  
  









};
