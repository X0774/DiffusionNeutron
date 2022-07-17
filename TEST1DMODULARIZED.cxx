#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;

//DEFINITION OF VARIABLES ON ALL SCOPES
double L,A,lambda,mu,eta;
int N;
int p=0;
double *a;

double Function (double x)
{
	return A*exp( -lambda*(pow(((x-L/2)/(L/2)),2) ));
};

double Int (double x)
{
	return 2*Function(x)*sin(p*M_PI*x/L)/L;
}

double NeutronDensity(double t,double x,double y,double z)
{
	double n=0;
	for (p=1;p<N;++p)
	{
		n=n+a[p]* exp ( eta *t - mu *pow(( p * M_PI / L ),2)* t ) * sin ( p * M_PI * x / L );
	};
	return n;
};

void plot(double (*func)(double, double, double, double), double tin, double tfin, double xin, double xfin, double yin, double yfin, double zin, double zfin, double dt, double dx)
{
	fstream newfile;
	int NPointsx=max(abs((xfin-xin)/dx),1.0); 
	int NPointsy=max(abs((yfin-yin)/dx),1.0);
	int NPointsz=max(abs((zfin-zin)/dx),1.0);
	int NStep=max(abs((tfin-tin)/dt),1.0);
	cout << tin<<"	"<<tfin<<"	"<<NStep;
	
	double n;
	newfile.open("output.txt",ios::out);
	newfile << "n(t,x,y,z)" <<"	"<<"t"<< "	" << "x" << "	" << "y" << "	" << "z" <<endl;
	for (int c=0; c<NStep; ++c)
		for (int i=0; i<NPointsx; ++i)
			for (int k=0; k<NPointsy; ++k)
				for (int l=0; l<NPointsz; ++l)
				{
					n=func(tin+c*dt,xin+i*dx,yin+k*dx,zin+l*dx);
					newfile << n << "	" << tin+c*dt<<"	"<< xin+i*dx << "	" << yin+k*dx << "	" << zin+l*dx << endl;
				};	
	newfile.close();
};

int main(int argc, char **argv)
{
	//Reading Parameters from a Configuration File
	fstream newfile;
	newfile.open("config.txt",ios::in);
	if (newfile.is_open())
	{
		newfile >> L >> A >> lambda >> N >> mu >> eta; //0.111 1 100 31 2.3446e+05 1.8958e+8 are the suggested values  
  	} else
  	{
   		exit(EXIT_FAILURE);
	}
	newfile.close(); //close the file object.
	
	
	//COMPUTE a[p]
	a=new double[N];
	Numerics num;
	for (p=0;p<N;++p)
	{
		a[p]=num.trapezoidal(Int,0,L,0.0001); 
		cout << a[p] << "	";
	}
	
	//Computing ns
	cout << "This is now a small test for x and t; please write dx, final t and dt" << endl;
	double dx; //0.001
	double tfin; //0.00002;
	double dt; //0.0000001;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> tfin;
	cout << "dt: ";
	cin >> dt;
	newfile.open("output.txt",ios::out);
	newfile << "n(x,t)" << "	" << "t" << "	" << "x" << endl;
	plot(NeutronDensity,0,tfin,0,L,0,0,0,0,dt,dx);
	newfile.close();
	
	delete[] a;
	return 0;
}
