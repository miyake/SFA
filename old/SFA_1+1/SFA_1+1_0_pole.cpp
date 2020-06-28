// SFA (1+1)サイト

#include <iostream.h>
#include <fstream.h>
#include <math.h>
//#include <process.h>
#include "dsyev.h"
#include <complex>
 typedef std::complex<double> Complex;


main()
{
	
	int i,j,m,n,nV,n_hh,n_ea,n_omega; int Nur,Ndr,Nug,Ndg; double Va,hha,eaa,Pa;
	int x,xx,k,kk; double K,KK;
	const double pi=3.1415926535897932; const double delta=0.01;
	// 生成消滅演算子
	double c1u_0[6][4],c1u_d1[4][6],c2u_0[6][4],c2u_d1[4][6];
	double c1d_0[6][4],c1d_u1[4][6],c2d_0[6][4],c2d_u1[4][6];
	for(i=0;i<=5;++i){
		for(j=0;j<=3;++j){
			c1u_0[i][j]=c1d_0[i][j]=c2u_0[i][j]=c2d_0[i][j]=c1u_d1[j][i]=c1d_u1[j][i]=c2u_d1[j][i]=c2d_u1[j][i]=0.0;
		};
	};
	 c1u_0[0][0]=c1u_0[3][2]=c1u_0[4][3] =c1u_d1[0][1]=c1u_d1[1][2]=c1u_d1[3][5]
	=c2u_0[0][1]=c2u_0[1][2]=c2u_0[2][3] =c2u_d1[0][3]=c2u_d1[1][4]=c2u_d1[2][5]
	=c1d_0[0][0]=c1d_0[2][2]=c1d_0[4][3] =c1d_u1[0][1]=c1d_u1[1][3]=c1d_u1[3][5]
	=c2d_0[0][1]=c2d_0[1][2]=c2d_0[3][3] =c2d_u1[0][2]=c2d_u1[1][4]=c2d_u1[2][5]
	=1.0;
	 c2u_0[2][3]=c2u_d1[0][3]=c1d_0[2][2]=c1d_u1[0][1]=c1d_u1[3][5]=c2d_u1[0][2]=c2d_u1[1][4]=c2d_u1[2][5]=-1.0;
	// 準備
	double e,t,U,h,mu,ec,ea,V,hh,P;
	double H0[6][6],Hu1[4][4],Hd1[4][4]; double Eu2,Ed2;
	double E0[6],Eu1[4],Ed1[4]; double work[191]; char jobz='V'; char uplo='U'; int info;
	double omega,Omega_r[25],Omega_t[25],Omega,Omega_a[25],Sur[25],Sdr[25],Sug[25],Sdg[25],Sr[25],Sg[25],Vaa[25];
	double gIRu1,gIRu2,gIRd1,gIRd2;
	double GRu11,GRu12,GRu21,GRu22,GRd11,GRd12,GRd21,GRd22,SEu11,SEu12,SEu21,SEu22,SEd11,SEd12,SEd21,SEd22;
	double Au11[5],Au12[5],Au21[5],Au22[5],Ad11[5],Ad12[5],Ad21[5],Ad22[5],Bu11[5],Bu12[5],Bu21[5],Bu22[5],Bd11[5],Bd12[5],Bd21[5],Bd22[5];
	double a,au1,au2,ad1,ad2,bu1,bu2,bd1,bd2;
	double gIu[64],gId[64];
	double nu[64],nd[64],Nu,Nd,N;
	
	Omega=0.0;
	// 元の系のパラメーター
	t=1.0;
	cout << "U="; cin >> U;
	cout << "mu="; cin >> mu;
	h=0.0;
	// 代理系のパラメーター
	ec=0.0;hh=0.0; P=0.0;
	cout << "V="; cin >> V;
	cout << "ea="; cin >> ea;
	if(V==0.0) V=0.0000001;
	Omega_t=Omega_r=Sur=Sdr=Sug=Sdg=0.0; Nur=Ndr=Nug=Ndg=0;
	for(i=0;i<=5;++i){
		for(j=0;j<=5;++j){
			H0[i][j]=0.0;
		};
	};
	for(i=0;i<=3;++i){
		for(j=0;j<=3;++j){
			Hu1[i][j]=Hd1[i][j]=0.0;
		};
	};
	// 代理系のハミルトニアンの対角成分
	H0[0][0]=0.0;
	Hu1[0][0]=(ec-mu)-hh; Hu1[1][1]=(ea-mu);
	Hd1[0][0]=(ec-mu)+hh; Hd1[1][1]=(ea-mu);
	H0[1][1]=(2*ec-2*mu)+U; H0[2][2]=(ec+ea-2*mu)-hh; H0[3][3]=(ec+ea-2*mu)+hh; H0[4][4]=(2*ea-2*mu);
	Eu2=(ec+ea-2*mu)-hh; Ed2=(ec+ea-2*mu)+hh;
	Hu1[2][2]=(2*ec+ea-3*mu)+U; Hu1[3][3]=(ec+2*ea-3*mu)-hh;
	Hd1[2][2]=(2*ec+ea-3*mu)+U; Hd1[3][3]=(ec+2*ea-3*mu)+hh;
	H0[5][5]=(2*ec+2*ea-4*mu)+U;
	// 代理系のホッピング
	Hu1[0][1]=Hu1[1][0] =Hd1[0][1]=Hd1[1][0] =H0[1][2]=H0[2][1]=H0[4][2]=H0[2][4] =V;
	H0[1][3]=H0[3][1]=H0[4][3]=H0[3][4] =Hu1[2][3]=Hu1[3][2] =Hd1[2][3]=Hd1[3][2] =-V;
	// 代理系のペアリング
	H0[0][1]=Hu1[1][2]=Hd1[1][2]=H0[4][5]
	=H0[1][0]=Hu1[2][1]=Hd1[2][1]=H0[5][4]
	=P;
	// 代理系のハミルトニアンの対角化
	dsyev_( jobz, uplo, 6, (double**)H0, 6, E0, work, 17, info);
	if(info!=0) cout << "info=" << info << " (H0)\n";
	dsyev_( jobz, uplo, 4, (double**)Hu1, 4, Eu1, work, 11, info);
	if(info!=0) cout << "info=" << info << " (Hu1)\n";
	dsyev_( jobz, uplo, 4, (double**)Hd1, 4, Ed1, work, 11, info);
	if(info!=0) cout << "info=" << info << " (Hd1)\n";
	if(E0[0]==E0[1]) cout << " H0 の基底状態は縮退しています\n";
	// 代理系の熱力学ポテンシャル
	Omega_r[n_ea]=E0[0]; //cout << "Omega_r=" << Omega_r << "\n";
	// 生成消滅演算子の基底のとりかえ、代理系のグリーン関数の重み、Rr
	//Rr=0.0;
	for(m=0;m<=3;++m){
		au1=bu1=au2=bu2=ad1=bd1=ad2=bd2=0.0;
		for(i=0;i<=5;++i){
			for(j=0;j<=3;++j){
				au1+=H0[0][i]*c1u_0[i][j]*Hu1[m][j]; bu1+=Hd1[m][j]*c1u_d1[j][i]*H0[0][i];
				au2+=H0[0][i]*c2u_0[i][j]*Hu1[m][j]; bu2+=Hd1[m][j]*c2u_d1[j][i]*H0[0][i];
				ad1+=H0[0][i]*c1d_0[i][j]*Hd1[m][j]; bd1+=Hu1[m][j]*c1d_u1[j][i]*H0[0][i];
				ad2+=H0[0][i]*c2d_0[i][j]*Hd1[m][j]; bd2+=Hu1[m][j]*c2d_u1[j][i]*H0[0][i];
			};
		};
		Au11[m]=au1*au1; Au12[m]=au1*au2; Au22[m]=au2*au2;
		Bu11[m]=bu1*bu1; Bu12[m]=bu1*bu2; Bu22[m]=bu2*bu2;
		Ad11[m]=ad1*ad1; Ad12[m]=ad1*ad2; Ad22[m]=ad2*ad2;
		Bd11[m]=bd1*bd1; Bd12[m]=bd1*bd2; Bd22[m]=bd2*bd2;
	};
	for(n_omega=0;n_omega<=199999;++n_omega){
	 omega=-20.0+0.0001*n_omega;
	 if((U==0.0)&&(fabs(omega-V)<0,0000000005)) omega+=0.000000001;
	 // 代理系のグリーン関数
	 GRu11=GRu12=GRu21=GRu22=GRd11=GRd12=GRd21=GRd22=0.0;
	 for(m=0;m<=3;++m){
	 	if(fabs(Au11[m])>1.0e-15) GRu11+=Au11[m]/(omega-Eu1[m]+E0[0]);
	 	if(fabs(Bu11[m])>1.0e-15) GRu11+=Bu11[m]/(omega-E0[0]+Ed1[m]);
	 	if(fabs(Au12[m])>1.0e-15) GRu12+=Au12[m]/(omega-Eu1[m]+E0[0]);
	 	if(fabs(Bu12[m])>1.0e-15) GRu12+=Bu12[m]/(omega-E0[0]+Ed1[m]);
	 	if(fabs(Au22[m])>1.0e-15) GRu22+=Au22[m]/(omega-Eu1[m]+E0[0]);
	 	if(fabs(Bu22[m])>1.0e-15) GRu22+=Bu22[m]/(omega-E0[0]+Ed1[m]);
	 	if(fabs(Ad11[m])>1.0e-15) GRd11+=Ad11[m]/(omega-Ed1[m]+E0[0]);
	 	if(fabs(Bd11[m])>1.0e-15) GRd11+=Bd11[m]/(omega-E0[0]+Eu1[m]);
	 	if(fabs(Ad12[m])>1.0e-15) GRd12+=Ad12[m]/(omega-Ed1[m]+E0[0]);
	 	if(fabs(Bd12[m])>1.0e-15) GRd12+=Bd12[m]/(omega-E0[0]+Eu1[m]);
	 	if(fabs(Ad22[m])>1.0e-15) GRd22+=Ad22[m]/(omega-Ed1[m]+E0[0]);
	 	if(fabs(Bd22[m])>1.0e-15) GRd22+=Bd22[m]/(omega-E0[0]+Eu1[m]);
	 };
	 // Rr
	 gIRu1=2/((GRu11+GRu22)+sqrt(pow((GRu11-GRu22),2)+4*pow(GRu12,2)));
	 gIRu2=2/((GRu11+GRu22)-sqrt(pow((GRu11-GRu22),2)+4*pow(GRu12,2)));
	 gIRd1=2/((GRd11+GRd22)+sqrt(pow((GRd11-GRd22),2)+4*pow(GRd12,2)));
	 gIRd2=2/((GRd11+GRd22)-sqrt(pow((GRd11-GRu22),2)+4*pow(GRd12,2)));
	 
	 if(gIRu1>0.0){
	 	Sur-=0.0001; Nur+=1;
	 };
	 if(gIRu2>0.0){
	 	Sur-=0.0001; Nur+=1;
	 };
	 if(gIRd1>0.0){
	 	Sdr-=0.0001; Ndr+=1;
	 };
	 if(gIRd2>0.0){
	 	Sdr-=0.0001; Ndr+=1;
	 };
	 // 代理系の自己エネルギー
	 SEu11=SEu12=SEu22=SEd11=SEd12=SEd22=0.0;
	 SEu11=omega+mu-ec+hh-GRu22/(GRu11*GRu22-GRu12*GRu12);
	 SEu12=-V+GRu12/(GRu11*GRu22-GRu12*GRu12);
	 SEu22=omega+mu-ea-GRu11/(GRu11*GRu22-GRu12*GRu12);
	 SEd11=omega+mu-ec-hh-GRd22/(GRd11*GRd22-GRd12*GRd12);
	 SEd12=-V+GRd12/(GRd11*GRd22-GRd12*GRd12);
	 SEd22=omega+mu-ea-GRd11/(GRd11*GRd22-GRd12*GRd12);
	 if((U==0.0)&&(fabs(SEu11)>0.000001)){
	 	cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu11="<<SEu11<<"\n";
	 	SEu11=0.0;
	 };
	 if((U==0.0)&&(fabs(SEd11)>0.000001)){
	 	cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd11="<<SEd11<<"\n";
	 	SEd11=0.0;
	 };
	 if(fabs(SEu12)>0.000001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu12="<<SEu12<<"\n";
	 if(fabs(SEu22)>0.000001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu22="<<SEu22<<"\n";
	 if(fabs(SEd12)>0.000001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd12="<<SEd12<<"\n";
	 if(fabs(SEu22)>0.000001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd22="<<SEd22<<"\n";
	 if((hh==0.0)&&(P=0.0)){
	 	a=SEu11-SEd11;
	 	if(fabs(a)>0.000001){
	 		cout << " hh=0,P=0 にもかかわらず SEu11!=SEd11 です。\n";
	 	};
	 	SEu11-=a/2.0; SEd11=SEu11;
	 };
	 // 元の系のグリーン関数の逆行列
	 for(k=0;k<=63;++k){
	 	K=pi*(k-31.0)/32.0;
	 	gIu[k]=omega+mu-e+h-SEu11+2*t*cos(K);
	 	gId[k]=omega+mu-e-h-SEd11+2*t*cos(K);
	 };
	 // Rg
	 for(k=0;k<=63;++k){
	 	if(gIu[k]>0.0){
	 		Sug[n_ea]-=0.0001/64.0; Nug+=1;
	 		//cout << "Nug+ : " << omega << "\n";
	 	};
	 	if(gId[k]>0.0){
	 		Sdg[n_ea]-=0.0001/64.0; Ndg+=1;
	 		//cout << "Ndg+ : " << omega << "\n";
	 	};
	 };
	 
	 
	 
	 
	}; // omega のループ、終わり
	Omega_t=Omega_r-Sur-Sdr+Sug+Sdg;
	Sr=Sur+Sdr; Sg=Sug+Sdg;
	cout <<"Omega_t="<<Omega_t<<"\n";
	cout <<"(Omega_r,Sr,Sg)=("<<Omega_r<<","<<Sr<<","<<Sg<<")\n";
	cout <<"(Sur,Sdr;Sug,Sdg)=("<<Sur<<","<<Sdg<<";"<<Sug<<","<<Sdg<<")\n";
	cout <<"(Nur,Ndr;Nug,Ndg)=("<<Nur<<","<<Ndr<<";"<<Nug<<","<<Ndg<<")\n";
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	cout << "\n";
	cout << "<< end! >>\n";
	
	
	
	
	
};
