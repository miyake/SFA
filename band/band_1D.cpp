// 1D, band

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <process.h>

main()
{
	int n,x,k; double K; const double pi=3.1415926535897932;
	double e,t,h; double band[64],bandu[64],bandd[64];
	
	t=1.0; e=0.0; h=0.1;
	
	for(k=0;k<=63;++k){
		K=pi*(k-31.0)/32.0;
		band[k]=e-2*t*cos(K);
		bandu[k]=e-h-2*t*cos(K);
		bandd[k]=e+h-2*t*cos(K);
	};
	
	ofstream fout;
	
	fout.open("band_1D.dat");
	for(k=0;k<=63;++k){
		fout << band[k] << "\n";
	};
	fout.close();
	
	fout.open("band_1D_u.dat");
	for(k=0;k<=63;++k){
		fout << bandu[k] << "\n";
	};
	fout.close();
	
	fout.open("band_1D_d.dat");
	for(k=0;k<=63;++k){
		fout << bandd[k] << "\n";
	};
	fout.close();
	
	
	
	
	
	
	
	
};
