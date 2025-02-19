// SFA, 2x2サイト

#include <iostream.h>
#include <math.h>
#include <process.h>

main()
{
	int i,j,k,n,N;
	double e,t,hF,hAF,U,mu;
	/*cout << "e="; cin >> e; cout << "t="; cin >> t; cout >> "U="; cin >> U;
	cout << "hF="; cin >> hF; cout << "hAF="; cin >> hAF; cout << "mu="; cin >> mu;*/
	e=0.0; t=10.0; hF=3.0; hAF=2.0; U=50.0; mu=20.0;
	// ハミルトニアンの用意
	double H00;
	double H10[4][4],H11[4][4];
	double H20[6][6],H21[16][16],H22[6][6];
	double H30[4][4],H31[24][24],H32[24][24],H33[4][4];
	double H40,H41[16][16],H42[36][36],H43[16][16],H44;
	double H51[4][4],H52[24][24],H53[24][24],H54[4][4];
	double H62[6][6],H63[16][16],H64[6][6];
	double H73[4][4],H74[4][4];
	double H84;
	for(i=0;i<=3;++i){
		for(j=0;j<=4;++j){
			H10[i][j]=0.0; H11[i][j]=0.0; H30[i][j]=0.0; H33[i][j]=0.0;
			H51[i][j]=0.0; H54[i][j]=0.0; H73[i][j]=0.0; H74[i][j]=0.0;
		};
	};
	for(i=0;i<=5;++i){
		for(j=0;j<=6;++j){
			H20[i][j]=0.0; H22[i][j]=0.0; H62[i][j]=0.0; H64[i][j]=0.0;
		};
	};
	for(i=0;i<=15;++i){
		for(j=0;j<=16;++j){
			H21[i][j]=0.0; H41[i][j]=0.0; H43[i][j]=0.0; H63[i][j]=0.0;
		};
	};
	for(i=0;i<=23;++i){
		for(j=0;j<=24;++j){
			H31[i][j]=0.0; H32[i][j]=0.0; H52[i][j]=0.0; H53[i][j]=0.0;
		};
	};
	for(i=0;i<=35;++i){
		for(j=0;j<=36;++j){
			H42[i][j]=0.0;
		};
	};
	// H00 & H40 & H44 & H84 (i=0)
	cout << "----- ( H00 & H40 & H44 & H84 ) --------\n";
	H00=0.0; H40=4*e-4*hF-4*mu; H44=4*e+4*hF-4*mu; H84=8*e+4*U-8*mu;
	cout <<"H00="<<H00<<", H40="<<H40<<", H44="<<H44<<", H84="<<H84<<"\n";
	// H10 & H11 (i=0,1,2,3)
	cout << "----------- ( H10 & H11 ) --------------\n";
	for(i=0;i<=3;++i){
		H10[i][i]=e-hF-mu; H11[i][i]=e+hF-mu;
	};
	H10[0][0]-=hAF; H10[1][1]+=hAF; H10[2][2]-=hAF; H10[3][3]+=hAF;
	H11[0][0]+=hAF; H11[1][1]-=hAF; H11[2][2]+=hAF; H11[3][3]-=hAF;
	H10[0][1]=H10[1][0]=H10[0][2]=H10[2][0]=H10[1][2]=H10[2][1]=H10[2][3]=H10[3][2]=t;
	H11[0][1]=H11[1][0]=H11[0][2]=H11[2][0]=H11[1][2]=H11[2][1]=H11[2][3]=H11[3][2]=t;
	for(i=0;i<=3;++i){
		cout <<"i="<<i<<", H10="<<H10[i][0]<<","<<H10[i][1]<<","<<H10[i][2]<<","<<H10[i][3]<<
		", // H11="<<H11[i][0]<<","<<H11[i][1]<<","<<H11[i][2]<<","<<H11[i][3]<<"\n";
	};
	// H20 & H22 (i=0,1,...,5)
	cout << "----------- ( H20 & H22 ) --------------\n";
	for(i=0;i<=5;++i){
		H20[i][i]=2*e-2*hF-2*mu; H22[i][i]=2*e+2*hF-2*mu;
	};
	H20[2][2]-=2*hAF; H20[3][3]+=2*hAF; H22[2][2]+=2*hAF; H22[3][3]-=2*hAF;
	H20[0][2]=H20[2][0]=H20[0][3]=H20[3][0]=H20[1][2]=H20[2][1]=H20[1][3]=H20[3][1]=t;
	H20[2][4]=H20[4][2]=H20[2][5]=H20[5][2]=H20[3][4]=H20[4][3]=H20[3][5]=H20[5][3]=t;
	H22[0][2]=H22[2][0]=H22[0][3]=H22[3][0]=H22[1][2]=H22[2][1]=H22[1][3]=H22[3][1]=t;
	H22[2][4]=H22[4][2]=H22[2][5]=H22[5][2]=H22[3][4]=H22[4][3]=H22[3][5]=H22[5][3]=t;
	for(i=0;i<=5;++i){
		cout <<"i="<<i<<", H20="<<H20[i][0]<<","<<H20[i][1]<<","<<H20[i][2]<<","<<H20[i][3]<<","<<H20[i][4]<<
		","<<H20[i][5]<<", // H22="<<H22[i][0]<<","<<H22[i][1]<<","<<H22[i][2]<<","<<H22[i][3]<<","<<H22[i][4]<<
		","<<H22[i][5]<<"\n";
	};
	// H21 (i=0,1,2,...,15)
	cout << "-------------- ( H21 ) -----------------\n";
	H21[0][0]=H21[5][5]=H21[10][10]=H21[15][15]=2*e+U-2*mu;
	H21[3][3]=H21[6][6]=H21[9][9]=H21[12][12]=2*e-2*mu;
	H21[1][1]=H21[2][2]=H21[13][13]=H21[14][14]=2*e-2*hAF-2*mu;
	H21[4][4]=H21[7][7]=H21[8][8]=H21[11][11]=2*e+2*hAF-2*mu;
	for(n=0;n<=3;++n){
		H21[4*n][4*n+1]=H21[4*n+1][4*n]=H21[4*n][4*n+2]=H21[4*n+2][4*n]=t;
		H21[4*n+1][4*n+3]=H21[4*n+3][4*n+1]=H21[4*n+2][4*n+3]=H21[4*n+3][4*n+2]=t;
	};
	for(i=0;i<=3;++i){
		H21[i][i+4]=H21[i+4][i]=H21[i][i+8]=H21[i+8][i]=t;
		H21[i+4][i+12]=H21[i+12][i+4]=H21[i+8][i+12]=H21[i+12][i+8]=t;
	};
	for(i=0;i<=15;++i){
		cout <<"i="<<i<<", H21="<<H21[i][0]<<","<<H21[i][1]<<","<<H21[i][2]<<","<<H21[i][3]<<", / "<<
		H21[i][4]<<","<<H21[i][5]<<","<<H21[i][6]<<","<<H21[i][7]<<", / "<<H21[i][8]<<","<<H21[i][9]<<
		","<<H21[i][10]<<","<<H21[i][11]<<", / "<<H21[i][12]<<","<<H21[i][13]<<","<<H21[i][14]<<","<<H21[i][15]<<"\n";
	};
	// H30 & H33 (i=0,1,2,3)
	cout << "----------- ( H30 & H33 ) --------------\n";
	for(i=0;i<=3;++i){
		H30[i][i]=3*e-3*hF-3*mu; H33[i][i]=3*e+3*hF-3*mu;
	};
	H30[0][0]+=hAF; H30[1][1]-=hAF; H30[2][2]-=hAF; H30[3][3]+=hAF;
	H33[0][0]-=hAF; H33[1][1]+=hAF; H33[2][2]+=hAF; H33[3][3]-=hAF;
	H30[0][1]=H30[1][0]=H30[0][2]=H30[2][0]=H30[1][3]=H30[3][1]=H30[2][3]=H30[3][2]=t;
	H33[0][1]=H33[1][0]=H33[0][2]=H33[2][0]=H33[1][3]=H33[3][1]=H33[2][3]=H33[3][2]=t;
	for(i=0;i<=3;++i){
		cout <<"i="<<i<<", H30="<<H30[i][0]<<","<<H30[i][1]<<","<<H30[i][2]<<","<<H30[i][3]<<
		", // H33="<<H33[i][0]<<","<<H33[i][1]<<","<<H33[i][2]<<","<<H33[i][3]<<"\n";
	};
	// H31 & H32 (i=0,1,2,...,23)
	cout << "----------- ( H31 & H32 ) --------------\n";
	for(i=0;i<=23;++i){
		H31[i][i]=3*e-hF-3*mu; H32[i][i]=3*e+hF-3*mu;
	};
	H31[0][0]+=U+hAF;    H31[1][1]+=U-hAF;   H31[2][2]+=-hAF;     H31[3][3]+=hAF;
	H31[4][4]+=U+hAF;    H31[5][5]+=-hAF;    H31[6][6]+=U-hAF;    H31[7][7]+=hAF;
	H31[8][8]+=U-hAF;    H31[9][9]+=-3*hAF;  H31[10][10]+=-3*hAF; H31[11][11]+=U-hAF;
	H31[12][12]+=3*hAF;  H31[13][13]+=U+hAF; H31[14][14]+=U+hAF;  H31[15][15]+=3*hAF;
	H31[16][16]+=hAF;    H31[17][17]+=U-hAF; H31[18][18]+=-hAF;   H31[19][19]+=U+hAF;
	H31[20][20]+=hAF;    H31[21][21]+=-hAF;  H31[22][22]+=U-hAF;  H31[23][23]+=U+hAF;
	H32[0][0]+=U-hAF;    H32[1][1]+=U+hAF;   H32[2][2]+=hAF;      H32[3][3]+=-hAF;
	H32[4][4]+=U-hAF;    H32[5][5]+=hAF;     H32[6][6]+=U+hAF;    H32[7][7]+=-hAF;
	H32[8][8]+=U+hAF;    H32[9][9]+=3*hAF;   H32[10][10]+=3*hAF;  H32[11][11]+=U+hAF;
	H32[12][12]+=-3*hAF; H32[13][13]+=U-hAF; H32[14][14]+=U-hAF;  H32[15][15]+=-3*hAF;
	H32[16][16]+=-hAF;   H32[17][17]+=U+hAF; H32[18][18]+=hAF;    H32[19][19]+=U-hAF;
	H32[20][20]+=-hAF;   H32[21][21]+=hAF;   H32[22][22]+=U+hAF;  H32[23][23]+=U-hAF;
	for(n=0;n<=5;++n){
		H31[4*n][4*n+1]=H31[4*n+1][4*n]=H31[4*n][4*n+2]=H31[4*n+2][4*n]=t;
		H31[4*n+1][4*n+3]=H31[4*n+3][4*n+1]=H31[4*n+2][4*n+3]=H31[4*n+3][4*n+2]=t;
		H32[4*n][4*n+1]=H32[4*n+1][4*n]=H32[4*n][4*n+2]=H32[4*n+2][4*n]=t;
		H32[4*n+1][4*n+3]=H32[4*n+3][4*n+1]=H32[4*n+2][4*n+3]=H32[4*n+3][4*n+2]=t;
	};
	for(i=0;i<=3;++i){
		H31[i][i+8]=H31[i+8][i]=H31[i][i+12]=H31[i+12][i]=H31[i+4][i+8]=H31[i+8][i+4]=H31[i+4][i+12]=H31[i+12][i+4]=t;
		H31[i+8][i+16]=H31[i+16][i+8]=H31[i+8][i+20]=H31[i+20][i+8]=H31[i+12][i+16]=H31[i+16][i+12]=H31[i+12][i+20]=H31[i+20][i+12]=t;
		H32[i][i+8]=H32[i+8][i]=H32[i][i+12]=H32[i+12][i]=H32[i+4][i+8]=H32[i+8][i+4]=H32[i+4][i+12]=H32[i+12][i+4]=t;
		H32[i+8][i+16]=H32[i+16][i+8]=H32[i+8][i+20]=H32[i+20][i+8]=H32[i+12][i+16]=H32[i+16][i+12]=H32[i+12][i+20]=H32[i+20][i+12]=t;
	};
	for(i=0;i<=23;++i){
		cout <<"i="<<i<<", H31="<<H31[i][0]<<","<<H31[i][1]<<","<<H31[i][2]<<","<<H31[i][3]<<",/"<<H31[i][4]<<","<<H31[i][5]<<","<<
		H31[i][6]<<","<<H31[i][7]<<",/"<<H31[i][8]<<","<<H31[i][9]<<","<<H31[i][10]<<","<<H31[i][11]<<",/"<<H31[i][12]<<","<<H31[i][13]<<
		","<<H31[i][14]<<","<<H31[i][15]<<",/"<<H31[i][16]<<","<<H31[i][17]<<","<<H31[i][18]<<","<<H31[i][19]<<",/"<<H31[i][20]<<","<<
		H31[i][21]<<","<<H31[i][22]<<","<<H31[i][23]<<"\n";
	};
	cout << "  - - - - - - - - - - - - - - - - - - - \n";
	for(i=0;i<=23;++i){
		cout <<"i="<<i<<", H32="<<H32[i][0]<<","<<H32[i][1]<<","<<H32[i][2]<<","<<H32[i][3]<<",/"<<H32[i][4]<<","<<H32[i][5]<<","<<
		H32[i][6]<<","<<H32[i][7]<<",/"<<H32[i][8]<<","<<H32[i][9]<<","<<H32[i][10]<<","<<H32[i][11]<<",/"<<H32[i][12]<<","<<H32[i][13]<<
		","<<H32[i][14]<<","<<H32[i][15]<<",/"<<H32[i][16]<<","<<H32[i][17]<<","<<H32[i][18]<<","<<H32[i][19]<<",/"<<H32[i][20]<<","<<
		H32[i][21]<<","<<H32[i][22]<<","<<H32[i][23]<<"\n";
	};
	// H41 & H43 (i=0,1,2,...,15)
	cout << "----------- ( H41 & H43 ) --------------\n";
	for(i=0;i<=15;++i){
		H41[i][i]=4*e-2*hF-4*mu; H43[i][i]=4*e+2*hF-4*mu;
	};
	H41[0][0]+=U+2*hAF; H41[1][1]+=U;       H41[2][2]+=U;         H41[3][3]+=2*hAF;
	H41[4][4]+=U;       H41[5][5]+=U-2*hAF; H41[6][6]+=-2*hAF;    H41[7][7]+=U;
	H41[8][8]+=U;       H41[9][9]+=-2*hAF;  H41[10][10]+=U-2*hAF; H41[11][11]+=U;
	H41[12][12]+=2*hAF; H41[13][13]+=U;     H41[14][14]+=U;       H41[15][15]+=U+2*hAF;
	H43[0][0]+=U+2*hAF; H43[1][1]+=U;       H43[2][2]+=U;         H43[3][3]+=2*hAF;
	H43[4][4]+=U;       H43[5][5]+=U-2*hAF; H43[6][6]+=-2*hAF;    H43[7][7]+=U;
	H43[8][8]+=U;       H43[9][9]+=-2*hAF;  H43[10][10]+=U-2*hAF; H43[11][11]+=U;
	H43[12][12]+=2*hAF; H43[13][13]+=U;     H43[14][14]+=U;       H43[15][15]+=U+2*hAF;
	for(n=0;n<=3;++n){
		H41[4*n][4*n+1]=H41[4*n+1][4*n]=H41[4*n][4*n+2]=H41[4*n+2][4*n]=t;
		H41[4*n+1][4*n+3]=H41[4*n+3][4*n+1]=H41[4*n+2][4*n+3]=H41[4*n+3][4*n+2]=t;
		H43[4*n][4*n+1]=H43[4*n+1][4*n]=H43[4*n][4*n+2]=H43[4*n+2][4*n]=t;
		H43[4*n+1][4*n+3]=H43[4*n+3][4*n+1]=H43[4*n+2][4*n+3]=H43[4*n+3][4*n+2]=t;
	};
	for(i=0;i<=3;++i){
		H41[i][i+4]=H41[i+4][i]=H41[i][i+8]=H41[i+8][i]=t;
		H41[i+4][i+12]=H41[i+12][i+4]=H41[i+8][i+12]=H41[i+12][i+8]=t;
		H43[i][i+4]=H43[i+4][i]=H43[i][i+8]=H43[i+8][i]=t;
		H43[i+4][i+12]=H43[i+12][i+4]=H43[i+8][i+12]=H43[i+12][i+8]=t;
	};
	for(i=0;i<=15;++i){
		cout <<"i="<<i<<", H41="<<H41[i][0]<<","<<H41[i][1]<<","<<H41[i][2]<<","<<H41[i][3]<<", / "<<
		H41[i][4]<<","<<H41[i][5]<<","<<H41[i][6]<<","<<H41[i][7]<<", / "<<H41[i][8]<<","<<H41[i][9]<<
		","<<H41[i][10]<<","<<H41[i][11]<<", / "<<H41[i][12]<<","<<H41[i][13]<<","<<H41[i][14]<<","<<H41[i][15]<<"\n";
	};
	cout << "  - - - - - - - - - - - - - - - - - - - \n";
	for(i=0;i<=15;++i){
		cout <<"i="<<i<<", H43="<<H43[i][0]<<","<<H43[i][1]<<","<<H43[i][2]<<","<<H43[i][3]<<", / "<<
		H43[i][4]<<","<<H43[i][5]<<","<<H43[i][6]<<","<<H43[i][7]<<", / "<<H43[i][8]<<","<<H43[i][9]<<
		","<<H43[i][10]<<","<<H43[i][11]<<", / "<<H43[i][12]<<","<<H43[i][13]<<","<<H43[i][14]<<","<<H43[i][15]<<"\n";
	};
	// H42 (i=0,1,2,...,35)
	cout << "-------------- ( H42 ) -----------------\n";
	for(i=0;i<=35;++i){
		H42[i][i]=4*e-4*mu;
	};
	H42[0][0]+=2*U;       H42[1][1]+=U;         H42[2][2]+=U+2*hAF;   H42[3][3]+=U-2*hAF;   H42[4][4]+=U;         H42[5][5]+=0.0;
	H42[6][6]+=U;         H42[7][7]+=2*U;       H42[8][8]+=U+2*hAF;   H42[9][9]+=U-2*hAF;   H42[10][10]+=0.0;     H42[11][11]+=U;
	H42[12][12]+=U-2*hAF; H42[13][13]+=U-2*hAF; H42[14][14]+=2*U;     H42[15][15]+=0.0;     H42[16][16]+=U-2*hAF; H42[17][17]+=U-2*hAF;
	H42[18][18]+=U+2*hAF; H42[19][19]+=U+2*hAF; H42[20][20]+=0.0;     H42[21][21]+=2*U;     H42[22][22]+=U+2*hAF; H42[23][23]+=U+2*hAF;
	H42[24][24]+=U;       H42[25][25]+=0.0;     H42[26][26]+=U+2*hAF; H42[27][27]+=U-2*hAF; H42[28][28]+=2*U;     H42[29][29]+=U;
	H42[30][30]+=0.0;     H42[31][31]+=U;       H42[32][32]+=U+2*hAF; H42[33][33]+=U-2*hAF; H42[34][34]+=U;       H42[35][35]+=2*U;
	for(n=0;n<=5;++n){
		H42[6*n][6*n+2]=H42[6*n+2][6*n]=H42[6*n][6*n+3]=H42[6*n+3][6*n]=t;
		H42[6*n+1][6*n+2]=H42[6*n+2][6*n+1]=H42[6*n+1][6*n+3]=H42[6*n+3][6*n+1]=t;
		H42[6*n+2][6*n+4]=H42[6*n+4][6*n+2]=H42[6*n+4][6*n+5]=H42[6*n+5][6*n+2]=t;
		H42[6*n+3][6*n+4]=H42[6*n+3][6*n+4]=H42[6*n+3][6*n+5]=H42[6*n+3][6*n+5]=t;
	};
	for(i=0;i<=5;++i){
		H42[i][i+12]=H42[i+12][i]=H42[i][i+18]=H42[i+18][i]=t;          H42[i+6][i+12]=H42[i+12][i+6]=H42[i+6][i+18]=H42[i+18][i+6]=t;
		H42[i+12][i+24]=H42[i+24][i+12]=H42[i+12][i+30]=H42[i+30][i+12]=H42[i+18][i+24]=H42[i+24][i+18]=H42[i+18][i+30]=H42[i+30][i+18]=t;
	};
	for(i=0;i<=35;++i){
		cout <<"i="<<i<<", H42="<<H42[i][0]<<","<<H42[i][1]<<","<<H42[i][2]<<","<<H42[i][3]<<","<<H42[i][4]<<","<<H42[i][5]<<",/"<<
		H42[i][6]<<","<<H42[i][7]<<","<<H42[i][8]<<","<<H42[i][9]<<","<<H42[i][10]<<","<<H42[i][11]<<",/"<<
		H42[i][12]<<","<<H42[i][13]<<","<<H42[i][14]<<","<<H42[i][15]<<","<<H42[i][16]<<","<<H42[i][17]<<",/"<<
		H42[i][18]<<","<<H42[i][19]<<","<<H42[i][20]<<","<<H42[i][21]<<","<<H42[i][22]<<","<<H42[i][23]<<",/"<<
		H42[i][24]<<","<<H42[i][25]<<","<<H42[i][26]<<","<<H42[i][27]<<","<<H42[i][28]<<","<<H42[i][29]<<",/"<<
		H42[i][30]<<","<<H42[i][31]<<","<<H42[i][32]<<","<<H42[i][33]<<","<<H42[i][34]<<","<<H42[i][35]<<"\n";
	};
	// H51 & H54 (i=0,1,2,3,4)
	cout << "----------- ( H51 & H54 ) --------------\n";
	H51[0][0]=H51[3][3]=5*e+U-3*hF+hAF-5*mu; H51[1][1]=H51[2][2]=5*e+U-3*hF-hAF-5*mu;
	H54[0][0]=H54[3][3]=5*e+U+3*hF-hAF-5*mu; H54[1][1]=H54[2][2]=5*e+U+3*hF+hAF-5*mu;
	H51[0][1]=H51[1][0]=H51[0][2]=H51[2][0]=H51[1][3]=H51[3][1]=H51[2][3]=H51[3][2]
	=H54[0][1]=H54[1][0]=H54[0][2]=H54[2][0]=H54[1][3]=H54[3][1]=H54[2][3]=H54[3][2]=t;
	for(i=0;i<=3;++i){
		cout <<"i="<<i<<", H51="<<H51[i][0]<<","<<H51[i][1]<<","<<H51[i][2]<<","<<H51[i][3]<<
		", // H54="<<H54[i][0]<<","<<H54[i][1]<<","<<H54[i][2]<<","<<H54[i][3]<<"\n";
	};
	// H52 & H53 (i=0,1,2,...,23)
	cout << "----------- ( H52 & H53 ) --------------\n";
	H52[0][0]=H52[1][1]=H52[8][8]  =H52[14][14]=H52[22][22]=H52[23][23]=2*U+hAF+(5*e-hF-5*mu);
	H52[3][3]=H52[6][6]=H52[10][10]=H52[13][13]=H52[17][17]=H52[20][20]=2*U-hAF+(5*e-hF-5*mu);
	H52[4][4]=H52[5][5]  =H52[18][18]=H52[19][19]=U+hAF+(5*e-hF-5*mu);
	H52[7][7]=H52[11][11]=H52[12][12]=H52[16][16]=U-hAF+(5*e-hF-5*mu);
	H52[2][2]=H52[20][20]=U+3*hAF+(5*e-hF-5*mu);
	H52[9][9]=H52[15][15]=U-3*hAF+(5*e-hF-5*mu);
	H53[0][0]=H53[1][1]=H53[8][8]  =H53[14][14]=H53[22][22]=H53[23][23]=2*U-hAF+(5*e+hF-5*mu);
	H53[3][3]=H53[6][6]=H53[10][10]=H53[13][13]=H53[17][17]=H53[20][20]=2*U+hAF+(5*e+hF-5*mu);
	H53[4][4]=H53[5][5]  =H53[18][18]=H53[19][19]=U-hAF+(5*e+hF-5*mu);
	H53[7][7]=H53[11][11]=H53[12][12]=H53[16][16]=U+hAF+(5*e+hF-5*mu);
	H53[2][2]=H53[20][20]=U-3*hAF+(5*e+hF-5*mu);
	H53[9][9]=H53[15][15]=U+3*hAF+(5*e+hF-5*mu);
	for(n=0;n<=3;++n){
		H52[6*n][6*n+2]=H52[6*n+2][6*n]=H52[6*n][6*n+3]=H52[6*n+3][6*n]
		=H52[6*n+1][6*n+2]=H52[6*n+2][6*n+1]=H52[6*n+1][6*n+3]=H52[6*n+3][6*n+1]
		=H52[6*n+2][6*n+4]=H52[6*n+4][6*n+2]=H52[6*n+2][6*n+5]=H52[6*n+5][6*n+2]
		=H52[6*n+3][6*n+4]=H52[6*n+4][6*n+3]=H52[6*n+3][6*n+5]=H52[6*n+5][6*n+3]
		=H53[6*n][6*n+2]=H53[6*n+2][6*n]=H53[6*n][6*n+3]=H53[6*n+3][6*n]
		=H53[6*n+1][6*n+2]=H53[6*n+2][6*n+1]=H53[6*n+1][6*n+3]=H53[6*n+3][6*n+1]
		=H53[6*n+2][6*n+4]=H53[6*n+4][6*n+2]=H53[6*n+2][6*n+5]=H53[6*n+5][6*n+2]
		=H53[6*n+3][6*n+4]=H53[6*n+4][6*n+3]=H53[6*n+3][6*n+5]=H53[6*n+5][6*n+3]=t;
	};
	for(i=0;i<=5;++i){
		H52[i][i+6]=H52[i+6][i]=H52[i][i+12]=H52[i+12][i]=H52[i+6][i+18]=H52[i+18][i+6]=H52[i+12][i+18]=H52[i+18][i+12]
		=H53[i][i+6]=H53[i+6][i]=H53[i][i+12]=H53[i+12][i]=H53[i+6][i+18]=H53[i+18][i+6]=H53[i+12][i+18]=H53[i+18][i+12]=t;
	};
	for(i=0;i<=23;++i){
		cout <<"i="<<i<<", H52="<<H52[i][0]<<","<<H52[i][1]<<","<<H52[i][2]<<","<<H52[i][3]<<","<<H52[i][4]<<","<<H52[i][5]<<",//"<<
		H52[i][6]<<","<<H52[i][7]<<","<<H52[i][8]<<","<<H52[i][9]<<","<<H52[i][10]<<","<<H52[i][11]<<",//"<<H52[i][12]<<","<<H52[i][13]<<
		","<<H52[i][14]<<","<<H52[i][15]<<","<<H52[i][16]<<","<<H52[i][17]<<",//"<<H52[i][18]<<","<<H52[i][19]<<","<<H52[i][20]<<","<<
		H52[i][21]<<","<<H52[i][22]<<","<<H52[i][23]<<"\n";
	};
	cout << "  - - - - - - - - - - - - - - - - - - - \n";
	for(i=0;i<=23;++i){
		cout <<"i="<<i<<", H53="<<H53[i][0]<<","<<H53[i][1]<<","<<H53[i][2]<<","<<H53[i][3]<<","<<H53[i][4]<<","<<H53[i][5]<<",//"<<
		H53[i][6]<<","<<H53[i][7]<<","<<H53[i][8]<<","<<H53[i][9]<<","<<H53[i][10]<<","<<H53[i][11]<<",//"<<H53[i][12]<<","<<H53[i][13]<<
		","<<H53[i][14]<<","<<H53[i][15]<<","<<H53[i][16]<<","<<H53[i][17]<<",//"<<H53[i][18]<<","<<H53[i][19]<<","<<H53[i][20]<<","<<
		H53[i][21]<<","<<H53[i][22]<<","<<H53[i][23]<<"\n";
	};
	// H62 & H64 (i=0,1,2,3,4,5)
	cout << "----------- ( H62 & H64 ) --------------\n";
	H62[0][0]=H62[1][1]=H62[4][4]=H62[5][5]=(6*e+2*U-2*hF-6*mu);
	H62[2][2]=(6*e+2*U-2*hF-6*mu)+2*hAF; H62[3][3]=(6*e+2*U-2*hF-6*mu)-2*hAF;
	H62[0][2]=H62[2][0]=H62[0][3]=H62[3][0]=H62[1][2]=H62[2][1]=H62[1][3]=H62[3][1]
	=H62[2][4]=H62[4][2]=H62[2][5]=H62[5][2]=H62[3][4]=H62[4][3]=H62[3][5]=H62[5][3]=t;
	H64[0][0]=H64[1][1]=H64[4][4]=H64[5][5]=(6*e+2*U+2*hF-6*mu);
	H64[2][2]=(6*e+2*U+2*hF-6*mu)-2*hAF; H64[3][3]=(6*e+2*U+2*hF-6*mu)+2*hAF;
	H64[0][2]=H64[2][0]=H64[0][3]=H64[3][0]=H64[1][2]=H64[2][1]=H64[1][3]=H64[3][1]
	=H64[2][4]=H64[4][2]=H64[2][5]=H64[5][2]=H64[3][4]=H64[4][3]=H64[3][5]=H64[5][3]=t;
	for(i=0;i<=5;++i){
		cout <<"i="<<i<<", H62="<<H62[i][0]<<","<<H62[i][1]<<","<<H62[i][2]<<","<<H62[i][3]<<","<<H62[i][4]<<","<<H62[i][5]<<
		", // H64="<<H64[i][0]<<","<<H64[i][1]<<","<<H64[i][2]<<","<<H64[i][3]<<","<<H64[i][4]<<","<<H64[i][5]<<"\n";
	};
	// H63 (i=0,1,2,...,15)
	cout << "-------------- ( H63 ) -----------------\n";
	H63[0][0]=H63[5][5]=H63[10][10]=H63[15][15]=(6*e-6*mu)+3*U;
	H63[1][1]=H63[2][2]=H63[13][13]=H63[14][14]=(6*e-6*mu)+2*U;
	H63[3][3]=H63[6][6]=H63[9][9]=H63[12][12]=(6*e-6*mu)+2*U+2*hAF;
	H63[4][4]=H63[7][7]=H63[8][8]=H63[11][11]=(6*e-6*mu)+2*U-2*hAF;
	for(n=0;n<=3;++n){
		H63[4*n][4*n+1]=H63[4*n+1][4*n]=H63[4*n][4*n+2]=H63[4*n+2][4*n]
		=H63[4*n+1][4*n+3]=H63[4*n+3][4*n+1]=H63[4*n+2][4*n+3]=H63[4*n+3][4*n+2]=t;
	};
	for(i=0;i<=3;++i){
		H63[i][i+4]=H63[i+4][i]=H63[i][i+8]=H63[i+8][i]
		=H63[i+4][i+12]=H63[i+12][i+4]=H63[i+8][i+12]=H63[i+12][i+8]=t;
	};
	for(i=0;i<=15;++i){
		cout <<"i="<<i<<", H63="<<H63[i][0]<<","<<H63[i][1]<<","<<H63[i][2]<<","<<H63[i][3]<<", / "<<
		H63[i][4]<<","<<H63[i][5]<<","<<H63[i][6]<<","<<H63[i][7]<<", / "<<H63[i][8]<<","<<H63[i][9]<<
		","<<H63[i][10]<<","<<H63[i][11]<<", / "<<H63[i][12]<<","<<H63[i][13]<<","<<H63[i][14]<<","<<H63[i][15]<<"\n";
	};
	// H73 & H74 (i=0,1,2,3)
	cout << "----------- ( H73 & H74 ) --------------\n";
	H73[0][0]=H73[3][3]=H74[1][1]=H74[2][2]=(7*e+3*U-hF-7*mu)-hAF;
	H73[1][1]=H73[2][2]=H74[0][0]=H74[3][3]=(7*e+3*U-hF-7*mu)+hAF;
	H73[0][1]=H73[1][0]=H73[0][2]=H73[2][0]=H73[1][2]=H73[2][1]=H73[2][3]=H73[3][2]=t;
	H74[0][1]=H74[1][0]=H74[0][2]=H74[2][0]=H74[1][2]=H74[2][1]=H74[2][3]=H74[3][2]=t;
	for(i=0;i<=3;++i){
		cout <<"i="<<i<<", H73="<<H73[i][0]<<","<<H73[i][1]<<","<<H73[i][2]<<","<<H73[i][3]<<
		", // H74="<<H74[i][0]<<","<<H74[i][1]<<","<<H74[i][2]<<","<<H74[i][3]<<"\n";
	};
	cout << "----------------------------------------\n";
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
