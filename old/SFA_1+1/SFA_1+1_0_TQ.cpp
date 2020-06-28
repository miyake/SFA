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
	double omega,Omega_r[31],Omega_t[31],Omega,Sur[31],Sdr[31],Sug[31],Sdg[31],Sr[31],Sg[31],gIRu1,gIRu2,gIRd1,gIRd2;
	double GRu11,GRu12,GRu21,GRu22,GRd11,GRd12,GRd21,GRd22,SEu11,SEu12,SEu21,SEu22,SEd11,SEd12,SEd21,SEd22;
	double Au11[5],Au12[5],Au21[5],Au22[5],Ad11[5],Ad12[5],Ad21[5],Ad22[5],Bu11[5],Bu12[5],Bu21[5],Bu22[5],Bd11[5],Bd12[5],Bd21[5],Bd22[5];
	double a,au1,au2,ad1,ad2,bu1,bu2,bd1,bd2;
	double gIu[64],gId[64];
	Complex omega_c;
	Complex GRu11c,GRu12c,GRu21c,GRu22c,GRd11c,GRd12c,GRd21c,GRd22c,SEu11c,SEu12c,SEu21c,SEu22c,SEd11c,SEd12c,SEd21c,SEd22c;
	double GRu11r,GRu12r,GRu21r,GRu22r,GRd11r,GRd12r,GRd21r,GRd22r,SEu11r,SEu12r,SEu21r,SEu22r,SEd11r,SEd12r,SEd21r,SEd22r;
	double GRu11i,GRu12i,GRu21i,GRu22i,GRd11i,GRd12i,GRd21i,GRd22i,SEu11i,SEu12i,SEu21i,SEu22i,SEd11i,SEd12i,SEd21i,SEd22i;
	Complex gIRu1c,gIRu2c,gIRd1c,gIRd2c,gIuc[64],gIdc[64],nuc[64],ndc[64];
	double gui,gdi,nu[64],nd[64],Nu,Nd,N;
	
	Omega=0.0;
	// 元の系のパラメーター
	t=1.0;
	cout << "U="; cin >> U;
	cout << "mu="; cin >> mu;
	h=0.0;
	// 代理系のパラメーター
	//cout << "V="; cin >> V;
	ec=0.0;hh=0.0; P=0.0;
	// hh のループ
	/*for(n_hh=0;n_hh<=2;++n_hh){
	hh=-0.01+0.01*n_hh;
	if(fabs(hh)<0.00001) hh=0.0;
	cout << "=======(hh="<<hh<<")======\n";*/
	// ea のループ
	//for(n_ea=0;n_ea<=4;++n_ea){
	//ea=U/2.0-0.02+0.01*n_ea;
	ea=U/2.0;
	cout << "------(ea="<<ea<<")-------\n";
	// V のループ
	for(nV=0;nV<=30;++nV){
	V=0.1*nV-1.5;
	if(nV==15) V=0.0000001;
	cout << "-------(V="<<V<<")-------\n";
	Omega_t[nV]=Omega_r[nV]=Sur[nV]=Sdr[nV]=Sug[nV]=Sdg[nV]=0.0; Nur=Ndr=Nug=Ndg=0;
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
	/*for(n=0;n<=5;++n){
		cout <<"n="<<n<<": E0[n]="<<E0[n]<<", H0[n][i]=("<<H0[n][0]<<","<<H0[n][1]<<","<<H0[n][2]<<","<<
		H0[n][3]<<","<<H0[n][4]<<","<<H0[n][5]<<")\n";
	};
	for(n=0;n<=3;++n){
		cout <<"n="<<n<<": Eu1[n]="<<Eu1[n]<<", Hu1[n][i]=("<<Hu1[n][0]<<","<<Hu1[n][1]<<","<<Hu1[n][2]<<","<<Hu1[n][3]<<")\n";
	};
	for(n=0;n<=3;++n){
		cout <<"n="<<n<<": Ed1[n]="<<Ed1[n]<<", Hd1[n][i]=("<<Hd1[n][0]<<","<<Hd1[n][1]<<","<<Hd1[n][2]<<","<<Hd1[n][3]<<")\n";
	};*/
	// 代理系の熱力学ポテンシャル
	Omega_r[nV]=E0[0]; //cout << "Omega_r=" << Omega_r << "\n";
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
		//cout <<"m="<<m<<": Au11="<<Au11[m]<<", Bu11="<<Bu11[m]<<"; Au12="<<Au12[m]<<", Bu12="<<Bu12[m]<<
		//", Au22="<<Au22[m]<<", Bu22="<<Bu22[m]<<"\n";
		//cout <<"m="<<m<<": Ad11="<<Ad11[m]<<", Bd11="<<Bd11[m]<<"; Ad12="<<Ad12[m]<<", Bd12="<<Bd12[m]<<
		//", Ad22="<<Ad22[m]<<", Bd22="<<Bd22[m]<<"\n";
		// Rr
		/*if(fabs((Au11[m])>0.0000001)&&(Eu1[m]-E0[0]<0)) Rr+=Eu1[m]-E0[0];
		if(fabs((Bu11[m])>0.0000001)&&(E0[0]-Ed1[0]<0)) Rr+=E0[0]-Ed1[m];
		if(fabs((Au12[m])>0.0000001)&&(Eu1[m]-E0[0]<0)) Rr+=Eu1[m]-E0[0];
		if(fabs((Bu12[m])>0.0000001)&&(E0[0]-Ed1[0]<0)) Rr+=E0[0]-Ed1[m];
		if(fabs((Au22[m])>0.0000001)&&(Eu1[m]-E0[0]<0)) Rr+=Eu1[m]-E0[0];
		if(fabs((Bu22[m])>0.0000001)&&(E0[0]-Ed1[0]<0)) Rr+=E0[0]-Ed1[m];
		if(fabs((Ad11[m])>0.0000001)&&(Ed1[m]-E0[0]<0)) Rr+=Ed1[m]-E0[0];
		if(fabs((Bd11[m])>0.0000001)&&(E0[0]-Eu1[0]<0)) Rr+=E0[0]-Eu1[m];
		if(fabs((Ad12[m])>0.0000001)&&(Ed1[m]-E0[0]<0)) Rr+=Ed1[m]-E0[0];
		if(fabs((Bd12[m])>0.0000001)&&(E0[0]-Eu1[0]<0)) Rr+=E0[0]-Eu1[m];
		if(fabs((Ad22[m])>0.0000001)&&(Ed1[m]-E0[0]<0)) Rr+=Ed1[m]-E0[0];
		if(fabs((Bd22[m])>0.0000001)&&(E0[0]-Eu1[0]<0)) Rr+=E0[0]-Eu1[m];*/
	};
	//cout << "Rr=" << Rr << "\n";
	for(n_omega=0;n_omega<=999999;++n_omega){
	 omega=-100.0+0.0001*n_omega;
	 if((U==0.0)&&(fabs(omega-V)<0,0000000005)) omega+=0.000000001;
	 //cout <<"------ ( omega = "<<omega<<" ) -------\n";
	 
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
	 //if(fabs(GRu11)<1.0e-15) GRu11=0.0;
	 //if(fabs(GRu12)<1.0e-15) GRu12=0.0;
	 //if(fabs(GRu22)<1.0e-15) GRu22=0.0;
	 //if(fabs(GRd11)<1.0e-15) GRd11=0.0;
	 //if(fabs(GRd12)<1.0e-15) GRd12=0.0;
	 //if(fabs(GRd22)<1.0e-15) GRd22=0.0;
	 // Rr
	 gIRu1=2/((GRu11+GRu22)+sqrt(pow((GRu11-GRu22),2)+4*pow(GRu12,2)));
	 gIRu2=2/((GRu11+GRu22)-sqrt(pow((GRu11-GRu22),2)+4*pow(GRu12,2)));
	 gIRd1=2/((GRd11+GRd22)+sqrt(pow((GRd11-GRd22),2)+4*pow(GRd12,2)));
	 gIRd2=2/((GRd11+GRd22)-sqrt(pow((GRd11-GRu22),2)+4*pow(GRd12,2)));
	 
	 if(gIRu1>0.0){
	 	Sur[nV]-=0.0001; Nur+=1;
	 };
	 if(gIRu2>0.0){
	 	Sur[nV]-=0.0001; Nur+=1;
	 };
	 if(gIRd1>0.0){
	 	Sdr[nV]-=0.0001; Ndr+=1;
	 };
	 if(gIRd2>0.0){
	 	Sdr[nV]-=0.0001; Ndr+=1;
	 };
	 //cout <<"GRu11="<<GRu11<<", GRu12="<<GRu12<<", GRu22="<<GRu22<<"\n";
	 //cout <<"GRd11="<<GRd11<<", GRd12="<<GRd12<<", GRd22="<<GRd22<<"\n";
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
	 //cout <<"SEu11="<<SEu11<<", SEu12="<<SEu12<<", SEu22="<<SEu22<<"\n";
	 //cout <<"SEd11="<<SEd11<<", SEd12="<<SEd12<<", SEd22="<<SEd22<<"\n";
	 // 元の系のグリーン関数の逆行列
	 for(k=0;k<=15;++k){
	 	K=pi*(k-7.0)/8.0;
	 	gIu[k]=omega+mu-e+h-SEu11+2*t*cos(K);
	 	gId[k]=omega+mu-e-h-SEd11+2*t*cos(K);
	 };
	 /*if(n_omega=0){
	 	for(k=0;k<=63;++k){
	 		cout <<"omega="<<omega<<", k="<<k<<": (gIu,gId)=("<<gIu[k]<<","<<gId[k]<<")\n";
	 	};
	 };*/
	 // Rg
	 for(k=0;k<=15;++k){
	 	if(gIu[k]>0.0){
	 		Sug[nV]-=0.0001/16.0; Nug+=1;
	 		//cout << "Nug+ : " << omega << "\n";
	 	};
	 	if(gId[k]>0.0){
	 		Sdg[nV]-=0.0001/16.0; Ndg+=1;
	 		//cout << "Ndg+ : " << omega << "\n";
	 	};
	 };
	 
	 
	 
	 
	}; // omega のループ、終わり
	Omega_t[nV]=Omega_r[nV]-Sur[nV]-Sdr[nV]+Sug[nV]+Sdg[nV]; Sr[nV]=Sur[nV]+Sdr[nV]; Sg[nV]=Sug[nV]+Sdg[nV];
	cout <<"Omega_t="<<Omega_t[nV]<<"\n";
	cout <<"(Omega_r,Sr,Sg)=("<<Omega_r[nV]<<","<<Sr[nV]<<","<<Sg[nV]<<")\n";
	cout <<"(Sur,Sdr;Sug,Sdg)=("<<Sur[nV]<<","<<Sdr[nV]<<";"<<Sug[nV]<<","<<Sdg[nV]<<")\n";
	cout <<"(Nur,Ndr;Nug,Ndg)=("<<Nur<<","<<Ndr<<";"<<Nug<<","<<Ndg<<")\n";
	if(Omega>Omega_t[nV]){
		Omega=Omega_t[nV]; Va=V; hha=hh; eaa=ea;
	};
	
	}; // V のループ、終わり
	//}; // hh のループ、終わり
	//}; // ea のループ、終わり
	cout <<"---------------------------\n";
	cout << "Va=" << Va << ", eaa=" << eaa << ", hha=" << hha << ", Omega=" << Omega << "\n";
	// 書き出し
	ofstream fout;
	fout.open("SFA_1+1_0_Omega_t.dat");
	for(nV=0;nV<=30;++nV){
		fout << Omega_t[nV] << "\n";
	};
	fout.close();
	
	fout.open("SFA_1+1_0_Omega_r.dat");
	for(nV=0;nV<=30;++nV){
		fout << Omega_r[nV] << "\n";
	};
	fout.close();
	
	fout.open("SFA_1+1_0_Sr.dat");
	for(nV=0;nV<=30;++nV){
		fout << Sr[nV] << "\n";
	};
	fout.close();
	
	fout.open("SFA_1+1_0_Sg.dat");
	for(nV=0;nV<=30;++nV){
		fout << Sg[nV] << "\n";
	};
	fout.close();
	
	
	
	
	
	// 最適化した後
	/*cout << "======= 最適化した後 =======\n";
	// パラメーター
	V=Va;
	// 初期化
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
	for(k=0;k<=63;++k){
		nuc[k]=ndc[k]=Complex(0.0,0.0);
		nu[k]=nd[k]=0.0;
	};
	Nu=Nd=N=0.0;
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
	Omega_r=E0[0]; //cout << "Omega_r=" << Omega_r << "\n";
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
		//cout <<"m="<<m<<": Au11="<<Au11[m]<<", Bu11="<<Bu11[m]<<"; Au12="<<Au12[m]<<", Bu12="<<Bu12[m]<<
		//", Au22="<<Au22[m]<<", Bu22="<<Bu22[m]<<"\n";
		//cout <<"m="<<m<<": Ad11="<<Ad11[m]<<", Bd11="<<Bd11[m]<<"; Ad12="<<Ad12[m]<<", Bd12="<<Bd12[m]<<
		//", Ad22="<<Ad22[m]<<", Bd22="<<Bd22[m]<<"\n";
	};
	//cout << "Rr=" << Rr << "\n";
	// omega のループ、始まり
	//cout <<"omega のループ、始まり\n";
	for(n=0;n<=999999;++n){
	 omega=-100.0+0.0001*n;
	 //if((U==0.0)&&(fabs(omega-V)<0,0000000005)) omega+=0.000000001;
	 omega_c=Complex(omega,delta);
	 //cout <<"------ ( omega = "<<omega<<" ) -------\n";
	 //cout << omega_c <<"\n";
	 // 代理系のグリーン関数
	 GRu11c=GRu12c=GRu21c=GRu22c=GRd11c=GRd12c=GRd21c=GRd22c=Complex(0.0,0.0);
	 for(m=0;m<=3;++m){
	 	if(fabs(Au11[m])>1.0e-15) GRu11c+=Complex(Au11[m],0.0)/(omega_c-Complex(Eu1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bu11[m])>1.0e-15) GRu11c+=Complex(Bu11[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Ed1[m],0.0));
	 	if(fabs(Au12[m])>1.0e-15) GRu12c+=Complex(Au12[m],0.0)/(omega_c-Complex(Eu1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bu12[m])>1.0e-15) GRu12c+=Complex(Bu12[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Ed1[m],0.0));
	 	if(fabs(Au22[m])>1.0e-15) GRu22c+=Complex(Au22[m],0.0)/(omega_c-Complex(Eu1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bu22[m])>1.0e-15) GRu22c+=Complex(Bu22[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Ed1[m],0.0));
	 	if(fabs(Ad11[m])>1.0e-15) GRd11c+=Complex(Ad11[m],0.0)/(omega_c-Complex(Ed1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bd11[m])>1.0e-15) GRd11c+=Complex(Bd11[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Eu1[m],0.0));
	 	if(fabs(Ad12[m])>1.0e-15) GRd12c+=Complex(Ad12[m],0.0)/(omega_c-Complex(Ed1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bd12[m])>1.0e-15) GRd12c+=Complex(Bd12[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Eu1[m],0.0));
	 	if(fabs(Ad22[m])>1.0e-15) GRd22c+=Complex(Ad22[m],0.0)/(omega_c-Complex(Ed1[m],0.0)+Complex(E0[0],0.0));
	 	if(fabs(Bd22[m])>1.0e-15) GRd22c+=Complex(Bd22[m],0.0)/(omega_c-Complex(E0[0],0.0)+Complex(Eu1[m],0.0));
	 };
	 // 代理系の自己エネルギー
	 SEu11c=omega_c+Complex(mu,0.0)-Complex(ec,0.0)+Complex(hh,0.0)-GRu22c/(GRu11c*GRu22c-GRu12c*GRu12c);
	 SEu12c=-Complex(V,0.0)+GRu12c/(GRu11c*GRu22c-GRu12c*GRu12c);
	 SEu22c=omega_c+Complex(mu,0.0)-Complex(ea,0.0)-GRu11c/(GRu11c*GRu22c-GRu12c*GRu12c);
	 SEd11c=omega_c+Complex(mu,0.0)-Complex(ec,0.0)-Complex(hh,0.0)-GRd22c/(GRd11c*GRd22c-GRd12c*GRd12c);
	 SEd12c=-Complex(V,0.0)+GRd12c/(GRd11c*GRd22c-GRd12c*GRd12c);
	 SEd22c=omega_c+(mu,0.0)-(ea,0.0)-GRd11c/(GRd11c*GRd22c-GRd12c*GRd12c);
	 // 元の系のグリーン関数の逆行列
	 for(k=0;k<=63;++k){
	 	K=pi*k/32.0-pi*31.0/32.0;
	 	gIuc[k]=omega_c+Complex(mu,0.0)-Complex(e,0.0)+Complex(h,0.0)-SEu11c+Complex(2*t*cos(K),0.0);
	 	gIdc[k]=omega_c+Complex(mu,0.0)-Complex(e,0.0)-Complex(h,0.0)-SEd11c+Complex(2*t*cos(K),0.0);
	 };
	 // １電子の分布関数
	 for(k=0;k<=63;++k){
	 	gui=((gIuc[k]).imag())/(((gIuc[k]).real())*((gIuc[k]).real())+((gIuc[k]).imag())*((gIuc[k]).imag()));
	 	gdi=((gIdc[k]).imag())/(((gIdc[k]).real())*((gIdc[k]).real())+((gIdc[k]).imag())*((gIdc[k]).imag()));
	 	nu[k]+=0.0001*gui/(64.0*pi); nd[k]+=0.0001*gdi/(64.0*pi);
	 };
	};// omega のループ、終わり
	
	// １電子の分布関数
	for(k=0;k<=63;++k){
		Nu+=nu[k]; Nd+=nd[k];
	};
	N=Nu+Nd;
	cout << "N=" << N << "\n";
	cout<<"(N;Nu.Nd)=("<<N<<";"<<Nu<<","<<Nd<<")\n";
	cout << "\n";
	cout << "k=(0,...,15)\n";
	cout<<"nu[k]=("<<nu[0]<<","<<nu[1]<<","<<nu[2]<<","<<nu[3]<<","<<nu[4]<<","<<nu[5]<<","<<nu[6]<<","<<nu[7]<<
	","<<nu[8]<<","<<nu[9]<<","<<nu[10]<<","<<nu[11]<<","<<nu[12]<<","<<nu[13]<<","<<nu[14]<<","<<nu[15]<<")\n";
	cout<<"nd[k]=("<<nd[0]<<","<<nd[1]<<","<<nd[2]<<","<<nd[3]<<","<<nd[4]<<","<<nd[5]<<","<<nd[6]<<","<<nd[7]<<
	","<<nd[8]<<","<<nd[9]<<","<<nd[10]<<","<<nd[11]<<","<<nd[12]<<","<<nd[13]<<","<<nd[14]<<","<<nd[15]<<")\n";
	cout << "k=(16,...,31)\n";
	cout<<"nu[k]=("<<nu[16]<<","<<nu[17]<<","<<nu[18]<<","<<nu[19]<<","<<nu[20]<<","<<nu[21]<<","<<nu[22]<<","<<nu[23]<<
	","<<nu[24]<<","<<nu[25]<<","<<nu[26]<<","<<nu[27]<<","<<nu[28]<<","<<nu[29]<<","<<nu[30]<<","<<nu[31]<<")\n";
	cout<<"nd[k]=("<<nd[16]<<","<<nd[17]<<","<<nd[18]<<","<<nd[19]<<","<<nd[20]<<","<<nd[21]<<","<<nd[22]<<","<<nd[23]<<
	","<<nd[24]<<","<<nd[25]<<","<<nd[26]<<","<<nd[27]<<","<<nd[28]<<","<<nd[29]<<","<<nd[30]<<","<<nd[31]<<")\n";
	cout << "k=(32,...,47)\n";
	cout<<"nu[k]=("<<nu[32]<<","<<nu[33]<<","<<nu[34]<<","<<nu[35]<<","<<nu[36]<<","<<nu[37]<<","<<nu[38]<<","<<nu[39]<<
	","<<nu[40]<<","<<nu[41]<<","<<nu[42]<<","<<nu[43]<<","<<nu[44]<<","<<nu[45]<<","<<nu[46]<<","<<nu[47]<<")\n";
	cout<<"nd[k]=("<<nd[32]<<","<<nd[33]<<","<<nd[34]<<","<<nd[35]<<","<<nd[36]<<","<<nd[37]<<","<<nd[38]<<","<<nd[39]<<
	","<<nd[40]<<","<<nd[41]<<","<<nd[42]<<","<<nd[43]<<","<<nd[44]<<","<<nd[45]<<","<<nd[46]<<","<<nd[47]<<")\n";
	cout << "k=(48,...,63)\n";
	cout<<"nu[k]=("<<nu[48]<<","<<nu[49]<<","<<nu[50]<<","<<nu[51]<<","<<nu[52]<<","<<nu[53]<<","<<nu[54]<<","<<nu[55]<<
	","<<nu[56]<<","<<nu[57]<<","<<nu[58]<<","<<nu[59]<<","<<nu[60]<<","<<nu[61]<<","<<nu[62]<<","<<nu[63]<<")\n";
	cout<<"nd[k]=("<<nd[48]<<","<<nd[49]<<","<<nd[50]<<","<<nd[51]<<","<<nd[52]<<","<<nd[53]<<","<<nd[54]<<","<<nd[55]<<
	","<<nd[56]<<","<<nd[57]<<","<<nd[58]<<","<<nd[59]<<","<<nd[60]<<","<<nd[61]<<","<<nd[62]<<","<<nd[63]<<")\n";
	*/
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	cout << "\n";
	cout << "<< end! >>\n";
	
	
	
	
	
};
