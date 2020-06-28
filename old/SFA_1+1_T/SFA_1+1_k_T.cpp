// SFA (1+1)サイト

#include <iostream.h>
#include <math.h>
#include <process.h>
#include "dsyev.h"

main()
{
	int i,j,m,n,N,n_omega,n_V,x,xx,k,kk; double K,KK; const double pi=3.141519;
	// 生成消滅演算子
	double c1u_0[6][4],c1u_u1[4][1],c1u_d1[4][6],c1u_d2[1][4];
	double c1d_0[6][4],c1d_u1[4][6],c1d_d1[4][1],c1d_u2[1][4];
	double c2u_0[6][4],c2u_u1[4][1],c2u_d1[4][6],c2u_d2[1][4];
	double c2d_0[6][4],c2d_u1[4][6],c2d_d1[4][1],c2d_u2[1][4];
	for(i=0;i<=5;++i){
		for(j=0;j<=3;++j){
			c1u_0[i][j]=c1d_0[i][j]=c2u_0[i][j]=c2d_0[i][j]
			=c1u_d1[j][i]=c1d_u1[j][i]=c2u_d1[j][i]=c2d_u1[j][i]
			=0.0;
		};
	};
	for(i=0;i<=3;++i){
		c1u_u1[i][0]=c1d_d1[i][0]=c2u_u1[i][0]=c2d_d1[i][0]
		=c1u_d2[0][i]=c1d_u2[0][i]=c2u_d2[0][i]=c2d_u2[0][i]
		=0.0;
	};
	 c1u_0[0][0]=c1u_0[3][2]=c1u_0[4][3] =c1u_u1[1][0]                           =c1u_d1[0][1]=c1u_d1[1][2]=c1u_d1[3][5] =c1u_d2[0][2]
	=c1d_0[0][0]=c1d_0[2][2]=c1d_0[4][3] =c1d_u1[0][1]=c1d_u1[1][3]=c1d_u1[3][5] =c1d_d1[1][0]                           =c1d_u2[0][2]
	=c2u_0[0][1]=c2u_0[1][2]=c2u_0[2][3] =c2u_u1[0][0]                           =c2u_d1[0][3]=c2u_d1[1][4]=c2u_d1[2][0] =c2u_d2[0][3]
	=c2d_0[0][1]=c2d_0[1][2]=c2d_0[3][3] =c2d_u1[0][2]=c2d_u1[1][4]=c2d_u1[2][5] =c2d_d1[0][0]                           =c2d_u2[0][3]
	=1.0;
	c1d_u1[0][1]=c1d_0[2][2]=c1d_u2[0][2]=c1d_u1[3][5]
	=c2u_u1[0][0]=c2u_d1[0][3]=c2u_0[2][3]=c2u_d2[0][3]
	=c2d_u1[0][2]=c2d_u1[1][4]=c2d_d1[0][0]=c2d_u1[2][5]
	=-1.0;
	// 準備
	double e,t,U,h,T,mu,ec,ea,V,hh,P; double Va;
	double H0[6][6],Hu1[4][4],Hd1[4][4]; double Eu2,Ed2;
	double E0[6],Eu1[4],Ed1[4]; double work[191]; char jobz='V'; char uplo='U'; int info;
	double C1[6][6],C2[6][6]; double a,Z,Omega_r,Omega_t,Omega,omega;
	double Au11,Au12,Au21,Au22,Ad11,Ad12,Ad21,Ad22;
	double GRu11,GRu12,GRu21,GRu22,GRd11,GRd12,GRd21,GRd22,SEu11,SEu12,SEu21,SEu22,SEd11,SEd12,SEd21,SEd22;
	double gIu[64],gId[64]; 
	double Rur,Rdr,Rug,Rdg; int Nur,Ndr,Nug,Ndg;
	// 元の系のパラメーター
	t=1.0;
	cout << "U="; cin >> U;
	cout << "T="; cin >> T;
	e=0.0; t=1.0; h=0.0; mu=U/2.0;
	// 代理系のパラメーター
	ec=0.0; ea=U/2.0; hh=0.0; P=0.0;
	// V のループ
	for(n_V=0;n_V<=12;++n_V){
	V=-1.2+0.2*n_V;
	cout << "-------(V="<<V<<")-------\n";
	Omega_t=Omega_r=Rur=Rdr=Rug=Rdg=0.0; Nur=Ndr=Nug=Ndg=0;
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
	for(i=0;i<=63;++i){
		gIu[i]=gId[i]=0.0;
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
	cout << "info=" << info << ", ";
	dsyev_( jobz, uplo, 4, (double**)Hu1, 4, Eu1, work, 11, info);
	cout << "info=" << info << ", ";
	dsyev_( jobz, uplo, 4, (double**)Hd1, 4, Ed1, work, 11, info);
	cout << "info=" << info << "\n";
	for(n=0;n<=5;++n){
		cout <<"n="<<n<<": E0[n]="<<E0[n]<<", H0[n][i]=("<<H0[n][0]<<","<<H0[n][1]<<","<<H0[n][2]<<","<<
		H0[n][3]<<","<<H0[n][4]<<","<<H0[n][5]<<")\n";
	};
	for(n=0;n<=3;++n){
		cout <<"n="<<n<<": Eu1[n]="<<Eu1[n]<<", Hu1[n][i]=("<<Hu1[n][0]<<","<<Hu1[n][1]<<","<<Hu1[n][2]<<","<<Hu1[n][3]<<")\n";
	};
	for(n=0;n<=3;++n){
		cout <<"n="<<n<<": Ed1[n]="<<Ed1[n]<<", Hd1[n][i]=("<<Hd1[n][0]<<","<<Hd1[n][1]<<","<<Hd1[n][2]<<","<<Hd1[n][3]<<")\n";
	};
	// 代理系の分配関数と熱力学ポテンシャル
	Z=0.0;
	for(n=0;n<=5;++n){
		Z+=exp(-E0[n]/T);
	};
	for(n=0;n<=3;++n){
		Z+=exp(-Eu1[n]/T)+exp(-Ed1[n]/T);
	};
	Omega_r=-T*log(Z);
	cout<<"Z="<<Z<<", Omega_r="<<Omega_r<<"\n";
	// 生成消滅演算子の表示のとりかえ
	for(n=0;n<=5;++n){
		for(m=0;m<=5;++m){
			C1[n][m]=C2[n][m]=0.0;
		};
	};
	// c1u_0, c2u_0
	for(n=0;n<=5;++n){
		for(m=0;m<=3;++m){
			for(i=0;i<=5;++i){
				for(j=0;j<=3;++j){
					C1[n][m]+=H0[n][i]*c1u_0[i][j]*Hu1[m][j]; C2[n][m]+=H0[n][i]*c2u_0[i][j]*Hu1[m][j]; 
				};
			};
			c1u_0[n][m]=C1[n][m]; c2u_0[n][m]=C2[n][m]; C1[n][m]=C2[n][m]=0.0;
		};
	};
	// c1u_u1, c2u_u2
	for(n=0;n<=3;++n){
			for(i=0;i<=3;++i){
					C1[n][0]+=Hu1[n][i]*c1u_u1[i][0]*Eu2; C2[n][0]+=Hu1[n][i]*c2u_u1[i][0]*Eu2;
			};
			c1u_u1[n][0]=C1[n][0]; c2u_u1[n][0]=C2[n][0]; C1[n][0]=C2[n][0]=0.0;
	};
	// c1u_d1, c2u_d1
	for(n=0;n<=3;++n){
		for(m=0;m<=5;++m){
			for(i=0;i<=3;++i){
				for(j=0;j<=5;++j){
					C1[n][m]+=Hd1[n][i]*c1u_d1[i][j]*H0[m][j]; C2[n][m]+=Hd1[n][i]*c2u_d1[i][j]*H0[m][j];
				};
			};
			c1u_d1[n][m]=C1[n][m]; c2u_d1[n][m]=C2[n][m]; C1[n][m]=C2[n][m]=0.0;
		};
	};
	// c1u_d2, c2u_d2
		for(m=0;m<=3;++m){
				for(j=0;j<=3;++j){
					C1[0][m]+=Ed2*c1u_d2[0][j]*Hd1[m][j]; C2[0][m]+=Ed2*c2u_d2[0][j]*Hd1[m][j];
				};
			c1u_d2[0][m]=C1[0][m]; c2u_d2[0][m]=C2[0][m]; C1[0][m]=C2[0][m]=0.0;
		};
	// c1d_0, c2d_0
	for(n=0;n<=5;++n){
		for(m=0;m<=3;++m){
			for(i=0;i<=5;++i){
				for(j=0;j<=3;++j){
					C1[n][m]+=H0[n][i]*c1d_0[i][j]*Hd1[m][j]; C2[n][m]+=H0[n][i]*c2d_0[i][j]*Hd1[m][j]; 
				};
			};
			c1d_0[n][m]=C1[n][m]; c2d_0[n][m]=C2[n][m]; C1[n][m]=C2[n][m]=0.0;
		};
	};
	// c1d_u1, c2d_u1
	for(n=0;n<=3;++n){
		for(m=0;m<=5;++m){
			for(i=0;i<=3;++i){
				for(j=0;j<=5;++j){
					C1[n][m]+=Hu1[n][i]*c1d_u1[i][j]*H0[m][j]; C2[n][m]+=Hu1[n][i]*c2d_u1[i][j]*H0[m][j]; 
				};
			};
			c1d_u1[n][m]=C1[n][m]; c2d_u1[n][m]=C2[n][m]; C1[n][m]=C2[n][m]=0.0;
		};
	};
	// c1d_d1, c2d_d1
	for(n=0;n<=3;++n){
			for(i=0;i<=3;++i){
					C1[n][0]+=Hd1[n][i]*c1d_d1[i][0]*Ed2; C2[n][m]+=Hd1[n][i]*c2d_d1[i][j]*Ed2;
			};
			c1d_d1[n][0]=C1[n][0]; c2d_d1[n][0]=C2[n][0]; C1[n][0]=C2[n][0]=0.0;
	};
	// c1d_u2, c2d_u2
		for(m=0;m<=3;++m){
				for(j=0;j<=3;++j){
					C1[0][m]+=Eu2*c1d_u2[0][j]*Hu1[m][j]; C2[0][m]+=Eu2*c2d_u2[0][j]*Hu1[m][j]; 
				};
			c1d_u2[0][m]=C1[0][m]; c2d_u2[0][m]=C2[0][m]; C1[0][m]=C2[0][m]=0.0;
		};
	// 代理系のグリーン関数
	for(n_omega=0;n_omega<=10;++n_omega){
	 omega=-0.5+0.1*n_omega;
	 GRu11=GRu12=GRu22=GRd11=GRd12=GRd22=0.0;
	 // u:(0,u1), d:(0,d1)
	 for(n=0;n<=5;++n){
	 	for(m=0;m<=3;++m){
	 		a=c1u_0[n][m]*c1u_0[n][m];
	 		if(fabs(a)>0.0000001) GRu11+=a*(exp(-E0[n]/T)+exp(-Eu1[m]/T))/(Z*(omega+E0[n]-Eu1[m]));
	 		a=c1u_0[n][m]*c2u_0[n][m];
	 		if(fabs(a)>0.0000001) GRu12+=a*(exp(-E0[n]/T)+exp(-Eu1[m]/T))/(Z*(omega+E0[n]-Eu1[m]));
	 		a=c2u_0[n][m]*c2u_0[n][m];
	 		if(fabs(a)>0.0000001) GRu22+=a*(exp(-E0[n]/T)+exp(-Eu1[m]/T))/(Z*(omega+E0[n]-Eu1[m]));
	 		a=c1d_0[n][m]*c1d_0[n][m];
	 		if(fabs(a)>0.0000001) GRd11+=a*(exp(-E0[n]/T)+exp(-Ed1[m]/T))/(Z*(omega+E0[n]-Ed1[m]));
	 		a=c1d_0[n][m]*c2d_0[n][m];
	 		if(fabs(a)>0.0000001) GRd12+=a*(exp(-E0[n]/T)+exp(-Ed1[m]/T))/(Z*(omega+E0[n]-Ed1[m]));
	 		a=c2d_0[n][m]*c2d_0[n][m];
	 		if(fabs(a)>0.0000001) GRd22+=a*(exp(-E0[n]/T)+exp(-Ed1[m]/T))/(Z*(omega+E0[n]-Ed1[m]));
	 	};
	 };
	 // u:(d1,0), d:(u1,0)
	 for(n=0;n<=3;++n){
	 	for(m=0;m<=5;++m){
	 		a=c1u_d1[n][m]*c1u_d1[n][m];
	 		if(fabs(a)>0.0000001) GRu11+=a*(exp(-Ed1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Ed1[n]-E0[m]));
	 		a=c1u_d1[n][m]*c2u_d1[n][m];
	 		if(fabs(a)>0.0000001) GRu12+=a*(exp(-Ed1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Ed1[n]-E0[m]));
	 		a=c2u_d1[n][m]*c2u_d1[n][m];
	 		if(fabs(a)>0.0000001) GRu22+=a*(exp(-Ed1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Ed1[n]-E0[m]));
	 		
	 		a=c1d_u1[n][m]*c1d_u1[n][m];
	 		if(fabs(a)>0.0000001) GRd11+=a*(exp(-Eu1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Eu1[n]-E0[m]));
	 		a=c1d_u1[n][m]*c2d_u1[n][m];
	 		if(fabs(a)>0.0000001) GRd12+=a*(exp(-Eu1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Eu1[n]-E0[m]));
	 		a=c2d_u1[n][m]*c2d_u1[n][m];
	 		if(fabs(a)>0.0000001) GRd22+=a*(exp(-Eu1[n]/T)+exp(-E0[m]/T))/(Z*(omega+Eu1[n]-E0[m]));
	 	};
	 };
	 // u:(u1,u2), d:(d1,d2)
	 for(n=0;n<=3;++n){
	 		a=c1u_u1[n][0]*c1u_u1[n][0];
	 		if(fabs(a)>0.0000001) GRu11+=a*(exp(-Eu1[n]/T)+exp(-Eu2/T))/(Z*(omega+Eu1[n]-Eu2));
	 		a=c1u_u1[n][0]*c2u_u1[n][0];
	 		if(fabs(a)>0.0000001) GRu12+=a*(exp(-Eu1[n]/T)+exp(-Eu2/T))/(Z*(omega+Eu1[n]-Eu2));
	 		a=c2u_u1[n][0]*c2u_u1[n][0];
	 		if(fabs(a)>0.0000001) GRu22+=a*(exp(-Eu1[n]/T)+exp(-Eu2/T))/(Z*(omega+Eu1[n]-Eu2));
	 		
	 		a=c1d_d1[n][0]*c1d_d1[n][0];
	 		if(fabs(a)>0.0000001) GRd11+=a*(exp(-Ed1[n]/T)+exp(-Ed2/T))/(Z*(omega+Ed1[n]-Ed2));
	 		a=c1d_d1[n][0]*c2d_d1[n][0];
	 		if(fabs(a)>0.0000001) GRd12+=a*(exp(-Ed1[n]/T)+exp(-Ed2/T))/(Z*(omega+Ed1[n]-Ed2));
	 		a=c2d_d1[n][0]*c2d_d1[n][0];
	 		if(fabs(a)>0.0000001) GRd22+=a*(exp(-Ed1[n]/T)+exp(-Ed2/T))/(Z*(omega+Ed1[n]-Ed2));
	 };
	 // u:(d2,d1), d:(u2,u1)
	 	for(m=0;m<=3;++m){
	 		a=c1u_d2[0][m]*c1u_d2[0][m];
	 		if(fabs(a)>0.0000001) GRu11+=a*(exp(-Ed2/T)+exp(-Ed1[m]))/(Z*(omega+Ed2-Ed1[m]));
	 		a=c1u_d2[0][m]*c2u_d2[0][m];
	 		if(fabs(a)>0.0000001) GRu12+=a*(exp(-Ed2/T)+exp(-Ed1[m]))/(Z*(omega+Ed2-Ed1[m]));
	 		a=c2u_d2[0][m]*c2u_d2[0][m];
	 		if(fabs(a)>0.0000001) GRu22+=a*(exp(-Ed2/T)+exp(-Ed1[m]))/(Z*(omega+Ed2-Ed1[m]));
	 		
	 		a=c1d_u2[0][m]*c1d_u2[0][m];
	 		if(fabs(a)>0.0000001) GRd11+=a*(exp(-Eu2/T)+exp(-Eu1[m]))/(Z*(omega+Eu2-Eu1[m]));
	 		a=c1d_u2[0][m]*c2d_u2[0][m];
	 		if(fabs(a)>0.0000001) GRd12+=a*(exp(-Eu2/T)+exp(-Eu1[m]))/(Z*(omega+Eu2-Eu1[m]));
	 		a=c2d_u2[0][m]*c2d_u2[0][m];
	 		if(fabs(a)>0.0000001) GRd22+=a*(exp(-Eu2/T)+exp(-Eu1[m]))/(Z*(omega+Eu2-Eu1[m]));
	 	};
	 //cout<<"GRu11="<<GRu11<<", GRu12="<<GRu12<<", GRu22="<<GRu22<<"\n";
	 //cout<<"GRd11="<<GRd11<<", GRd12="<<GRd12<<", GRd22="<<GRd22<<"\n";
	 // 代理系の自己エネルギー
	 SEu11=omega+mu-ec+hh-GRu22/(GRu11*GRu22-GRu12*GRu12);
	 SEu12=-V+GRu12/(GRu11*GRu22-GRu12*GRu12);
	 SEu22=omega+mu-ea-GRu11/(GRu11*GRu22-GRu12*GRu12);
	 SEd11=omega+mu-ec-hh-GRd22/(GRd11*GRd22-GRd12*GRd12);
	 SEd12=-V+GRd12/(GRd11*GRd22-GRd12*GRd12);
	 SEd22=omega+mu-ea-GRd11/(GRd11*GRd22-GRd12*GRd12);
	 //cout<<"SEu11="<<SEu11<<", SEu12="<<SEu12<<", SEu22="<<SEu22<<"\n";
	 //cout<<"SEd11="<<SEd11<<", SEd12="<<SEd12<<", SEd22="<<SEd22<<"\n";
	 if(fabs(SEu11)<0.00001) SEu11=0.0;
	 if(fabs(SEd11)<0.00001) SEd11=0.0;
	 if((U==0.0)&&(fabs(SEu11)>0.00001)) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu11="<<SEu11<<"\n";
	 if((U==0.0)&&(fabs(SEd11)>0.00001)) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd11="<<SEd11<<"\n";
	 if(fabs(SEu12)>0.00001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu12="<<SEu12<<"\n";
	 if(fabs(SEu22)>0.00001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEu22="<<SEu22<<"\n";
	 if(fabs(SEd12)>0.00001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd12="<<SEd12<<"\n";
	 if(fabs(SEu22)>0.00001) cout <<"(V,omega)=("<<V<<","<<omega<<") : SEd22="<<SEd22<<"\n";
	 if((hh==0.0)&&(P=0.0)){
	 	a=SEu11-SEd11;
	 	if(fabs(a)>0.00001){
	 		cout << " hh=0,P=0 にもかかわらず SEu11!=SEd11 です。\n";
	 	};
	 	SEu11-=a/2.0; SEd11=SEu11;
	 };
	 for(k=0;k<=63;++k){
	 	K=pi*k/32.0-pi*31.0/32.0;
	 	gIu[k]=omega+mu-e+h-SEu11+2*t*cos(K);
	 	gId[k]=omega+mu-e-h-SEd11+2*t*cos(K);
	 };
	 // Rg
	 for(k=0;k<=63;++k){
	 	if(fabs(gIu[k])<0.0005){
	 		Rug+=omega/128.0; Nug+=1;
	 		//cout << "Nug+ : " << omega << "\n";
	 	};
	 	if(fabs(gId[k])<0.0005){
	 		Rdg+=omega/128.0; Ndg+=1;
	 		//cout << "Ndg+ : " << omega << "\n";
	 	};
	 };
	 
	 
	 
	 
	}; // omega のループ、終わり
	//cout << "Rg=" << Rg << "\n";
	Omega_t=Omega_r-Rur-Rdr+Rug+Rdg;
	cout <<"Omega_t="<<Omega_t<<"\n";
	cout <<"(Omega_r/Rur.Rdr/Rug.Rdg)=("<<Omega_r<<"/"<<Rur<<","<<Rdr<<"/"<<Rug<<","<<Rdg<<")\n";
	cout <<"(Nur,Ndr;Nug,Ndg)=("<<Nur<<","<<Ndr<<";"<<Nug<<","<<Ndg<<")\n";
	if(Omega>Omega_t){
		Omega=Omega_t; Va=V;
	};
	
	}; // V のループ、終わり
	cout << "---------------------------\n";
	cout << "Va=" << Va << ", Omega=" << Omega << "\n";
	
	cout << "\n";
	cout << " << end! >>\n";
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
};
