// 1D, band

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <process.h>

main()
{
	int n,x,y,kx,ky; double Kx,Ky; const double pi=3.1415926535897932;
	double e,t,tn,h; double band[64][64],bandu[64][64],bandd[64][64];
	
	t=1.0; e=0.0; h=0.1;
	cout << "tn="; cin >> tn;
	
	for(kx=0;kx<=63;++kx){
		for(ky=0;ky<=63;++ky){
			Kx=pi*(kx-31.0)/32.0; Ky=pi*(ky-31.0)/32.0;
			band[kx][ky]=e-2*t*cos(Kx)-2*t*cos(Kx)-4*tn*cos(Kx)*cos(Ky);
		};
	};
	
	ofstream fout;
	
	fout.open("band_2D_kx.dat");
	for(kx=0;kx<=63;++kx){
		fout << band[kx][31] << "\n";
	};
	fout.close();
	
	fout.open("band_2D_X.dat");
	for(n=0;n<=32;++n){
		fout << band[63-n][31] << "\n";
	};
	for(n=0;n<=32;++n){
		fout << band[31+n][31+n] << "\n";
	};
	for(n=0;n<=32;++n){
		fout << band[63][63-n] << "\n";
	};
	for(n=0;n<=16;++n){
		fout << band[63-n][31+n] << "\n";
	};
	
	
	fout.close();
	
	
	
	
	
	
	
};
