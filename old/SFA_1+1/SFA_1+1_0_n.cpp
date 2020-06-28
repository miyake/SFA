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
	
	int i,j,m,n,nn,n_omega,Nur,Ndr,Nug,Ndg;
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
	double H0[6][6],Hu1[4][4],Hd1[4][4];
	double E0[6],Eu1[4],Ed1[4]; double work[191]; char jobz='V'; char uplo='U'; int info;
	double omega,Omega_r,Omega_t,Omega,Sur,Sdr,Sug,Sdg,Sr,Sg;
	double Au11[5],Au12[5],Au21[5],Au22[5],Ad11[5],Ad12[5],Ad21[5],Ad22[5],Bu11[5],Bu12[5],Bu21[5],Bu22[5],Bd11[5],Bd12[5],Bd21[5],Bd22[5];
	double a,au1,au2,ad1,ad2,bu1,bu2,bd1,bd2;
	Complex omega_c;
	Complex GRu11c,GRu12c,GRu21c,GRu22c,GRd11c,GRd12c,GRd21c,GRd22c,SEu11c,SEu12c,SEu21c,SEu22c,SEd11c,SEd12c,SEd21c,SEd22c;
	double GRu11r,GRu12r,GRu21r,GRu22r,GRd11r,GRd12r,GRd21r,GRd22r,SEu11r,SEu12r,SEu21r,SEu22r,SEd11r,SEd12r,SEd21r,SEd22r;
	double GRu11i,GRu12i,GRu21i,GRu22i,GRd11i,GRd12i,GRd21i,GRd22i,SEu11i,SEu12i,SEu21i,SEu22i,SEd11i,SEd12i,SEd21i,SEd22i;
	Complex gIRu1c,gIRu2c,gIRd1c,gIRd2c,gIuc[64],gIdc[64];
	double gui,gdi,nu[64],nd[64],Nu,Nd,N;
	// 元の系のパラメーター
	t=1.0;
	cout << "U="; cin >> U;
	mu=U/2.0;
	h=0.0; //cout << "h="; cin >> h;
	// 代理系のパラメーター
	cout << "V="; cin >> V;
	ea=U/2.0; //cout << "ea="; cin >> ea;
	hh=0.0; //cout << "hh="; cin >> hh;
	ec=0.0; P=0.0;
	if(V==0.0) V=0.0000001;
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
		nu[k]=nd[k]=0.0;
	};
	Nu=Nd=N=0.0;
	// 代理系のハミルトニアンの対角成分
	H0[0][0]=0.0;
	Hu1[0][0]=(ec-mu)-hh; Hu1[1][1]=(ea-mu);
	Hd1[0][0]=(ec-mu)+hh; Hd1[1][1]=(ea-mu);
	H0[1][1]=(2*ec-2*mu)+U; H0[2][2]=(ec+ea-2*mu)-hh; H0[3][3]=(ec+ea-2*mu)+hh; H0[4][4]=(2*ea-2*mu);
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
	// 生成消滅演算子の基底のとりかえ、代理系のグリーン関数の重み
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
	// omega のループ、始まり
	for(n=0;n<=199999;++n){
	 omega=-20.0+0.0001*n;
	 //if((U==0.0)&&(fabs(omega-V)<0,0000000005)) omega+=0.000000001;
	 omega_c=Complex(omega,delta);
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
	 // 元の系のグリーン関数の虚部、１電子の分布関数
	 for(k=0;k<=63;++k){
	 	gui=((gIuc[k]).imag())/(((gIuc[k]).real())*((gIuc[k]).real())+((gIuc[k]).imag())*((gIuc[k]).imag()));
	 	gdi=((gIdc[k]).imag())/(((gIdc[k]).real())*((gIdc[k]).real())+((gIdc[k]).imag())*((gIdc[k]).imag()));
	 	nu[k]+=0.0001*gui/pi; nd[k]+=0.0001*gdi/pi;
	 };
	};// omega のループ、終わり
	
	// 電子数
	for(k=0;k<=63;++k){
		Nu+=nu[k]/64.0; Nd+=nd[k]/64.0;
	};
	N=Nu+Nd;
	cout << "N=" << N << "\n";
	cout<<"(N;Nu.Nd)=("<<N<<";"<<Nu<<","<<Nd<<")\n";
	cout << "\n";
	
	ofstream fout;
	fout.open("SFA_1+1_0_nu.dat");
	for(k=0;k<=63;++k){
		fout << nu[k] << "\n";
	};
	fout.close();
	
	fout.open("SFA_1+1_0_nd.dat");
	for(k=0;k<=63;++k){
		fout << nd[k] << "\n";
	};
	fout.close();
	
	
	
	
	
	
	
	cout << "\n";
	cout << "<< end! >>\n";
	
	
	
	
	
};
