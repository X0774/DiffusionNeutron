#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;
//DEFINITION OF VARIABLES
double L,A,lambda,mu,eta;
int N;
int p=0;

double Function (double x)
{
	return A*exp( -lambda*(pow(((x-L/2)/(L/2)),2) ));
};

double Int (double x)
{
	return 2*Function(x)*sin(p*M_PI*x/L)/L;
}


int main(int argc, char **argv)
{
	//Reading Parameters from a Configuration File
	fstream newfile;
	newfile.open("config.txt",ios::in);
	if (newfile.is_open())
	{
	  newfile >> L >> A >> lambda >> N >> mu >> eta; //0.111 1 100 30 2.3446e+05 1.8958e+8 are the suggested values
      
    } else
    {
		exit(EXIT_FAILURE);
	}
      newfile.close(); //close the file object.
	
	
	
	double *a=new double[N];
	Numerics num;
	double n;
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
	double t,x;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> tfin;
	cout << "dt: ";
	cin >> dt;
	int NPoints=abs(L/dx);
	int NStep=abs(tfin/dt);
	newfile.open("output.txt",ios::out);
	newfile << "n(x,t)" << "	" << "t" << "	" << "x" << endl;
	for (int c=1; c<NStep; ++c)
	{
		t=c*dt;
		for (int i=1; i<NPoints; ++i)
		{
			x=i*dx;
			n=0; //At the end of every loop in x, n must be reset
			for (p=1;p<N;++p)
			{
				//n=n+a[p]*exp(  eta*c*dt - mu*pow((p*M_PI/L),2)*c*dt )*sin( p*M_PI*i*dx / L );
				n=n+a[p]* exp ( eta *t - mu *pow(( p * M_PI / L ),2)* t ) * sin ( p * M_PI * x / L );
			};
			newfile << n << "	" << t << "	" << 0+x << endl;	
		};
	};
	newfile.close();
	
	delete[] a;
	return 0;
}

