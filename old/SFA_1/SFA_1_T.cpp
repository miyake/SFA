// SFA, 1�T�C�g, T>0

#include <iostream.h>
#include <math.h>
#include <process.h>
#include "dsyev.h"

main()
{

int i,j,k,l,m,n,N,n_omega,Npr;
// �������ŉ��Z�q
double cu[4][4],cd[4][4];
for(i=0;i<=3;++i){
	for(j=0;j<=3;++j){
		cu[i][j]=cd[i][j]=0.0;
	};
};
cu[0][1]=cu[2][3]=cd[0][2]=1.0; cd[1][3]=-1.0;
// ����
double e,t,U,h, T,mu, ee,hh,P;
double H[4][4]; double E[4],work[191]; char jobz='V'; char uplo='U'; int info;
double Z,Omega_r,Rr,Rg,Omega_t,Omega, omega,a; double au[4][4],ad[4][4];
double GRu,GRd,SEu,SEd; double GIu[64][64],GId[64][64],gIu[64],gId[64];
// ���̌n�̃p�����[�^�[
e=0.0; t=1.0; h=0.0;
cout << "U="; cin >> U;
cout << "T="; cin >> T;
cout << "mu="; cin >> mu;
// �㗝�n�̃p�����[�^�[
ee=0.0; hh=0.0; P=0.0;
// ������
for(i=0;i<=3;++i){
	for(j=0;j<=3;++j){
		H[i][j]=au[i][j]=ad[i][j]=0.0;
	};
};
Z=Omega_r=Rr=Rg=0.0; Npr=0;
// �㗝�n�̃n�~���g�j�A��
H[0][0]=0.0; H[1][1]=ee-hh-mu; H[2][2]=ee+hh-mu; H[3][3]=2*ee+U-2*mu;
H[0][3]=H[3][0]=P;
// �㗝�n�̃n�~���g�j�A���̑Ίp��
dsyev_( jobz, uplo, 4, (double**)H, 4, E, work, 11, info);
if(info!=0) cout << "info=" << info << "\n";
// �㗝�n�̕��z�֐��ƔM�͊w�|�e���V����
for(i=0;i<=3;++i){
	Z+=exp(-E[i]/T);
};
Omega_r=-T*log(Z);
// �������ŉ��Z�q�̊��̂Ƃ肩���ARr
for(n=0;n<=3;++n){
	for(m=0;m<=3;++m){
		for(i=0;i<=3;++i){
			for(j=0;j<=3;++j){
				au[n][m]+=H[n][i]*cu[i][j]*H[m][j];
				ad[n][m]+=H[n][i]*cd[i][j]*H[m][j];
			};
		};
		if(fabs(au[n][m])<0.00001){
			Rr+=-T*log(1+exp((E[n]-E[m])/T)); Npr+=1;
		};
		if(fabs(ad[n][m])<0.00001){
			Rr+=-T*log(1+exp((E[n]-E[m])/T)); Npr+=1;
		};
	};
};
// omega �X�^�[�g
for(n_omega=0;n_omega<=500;++n_omega){
  omega=-5.0+0.01*n_omega;
  // �㗝�n�̃O���[���֐�
  GRu=GRd=0.0;
  for(n=0;n<=3;++n){
  	for(m=0;m<=3;++m){
  		if(fabs(au[n][m])>0.0000001) GRu+=au[n][m]*au[n][m]*(exp(-E[n]/T)+exp(-E[m]/T))/(omega+E[n]-E[m]);
  		if(fabs(ad[n][m])>0.0000001) GRd+=ad[n][m]*ad[n][m]*(exp(-E[n]/T)+exp(-E[m]/T))/(omega+E[n]-E[m]);
  	};
  };
  // �㗝�n�̎��ȃG�l���M�[
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
}; // omega �̃��[�v�A�I���
// Omega_t
Omega_t=Omega_r-Rr+Rg;
cout <<"Omega_t="<<Omega_t<<" ( Omega_r="<<Omega_r<<", Rr="<<Rr<<", Rg=" << Rg << " )\n";
  
  
  
  
  









};
