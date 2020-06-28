// report


#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <process.h>

main()
{
	double pi=3.1415,h=1.0545,c=2.9979,e=4.80325;
	double a_Au=4.08,a_Ag=4.09,a_Cu=3.61;
	double delta_Au,delta_Ag,delta_Cu;
	
	delta_Au=(pow(12*pi*pi,2/3)*h*c)/(2*e*a_Au*a_Au);
	delta_Ag=(pow(12*pi*pi,2/3)*h*c)/(2*e*a_Ag*a_Ag);
	delta_Cu=(pow(12*pi*pi,2/3)*h*c)/(2*e*a_Cu*a_Cu);
	
	cout << "Au : " << delta_Au << "\n";
	cout << "Ag : " << delta_Ag << "\n";
	cout << "Cu : " << delta_Cu << "\n";
	
	
	
};
