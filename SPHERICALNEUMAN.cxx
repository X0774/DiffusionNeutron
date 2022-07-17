#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"
using namespace std;

//DEFINING VARIABLES THAT MUST BE IN ALL SCOPES
double alpha,eta,mu,lambdat,r1,R;

double CriticalRadius(double r)
{
	return -1 + 3*r/lambdat/2 + r*sqrt(eta/mu)*cos(r*sqrt(eta/mu))/sin(r*sqrt(eta/mu));
};

double NewAlpha (double a)
{
	//return -1 + 3*(R/lambdat)/2 + R*sqrt((eta+a)/mu)*cos(R*sqrt((eta+a)/mu))/sin(R*sqrt((eta+a)/mu)); 
	return -1 + sqrt((eta+a)/mu) * R * (cos (R* sqrt((eta+a)/mu))/ sin ( R*sqrt((eta+a)/mu)) ) +(3/2) * R / lambdat;
};

double AlphaPrime(double r)
{
	Numerics num;
	return num.derivative(NewAlpha, 0.00001, r);
};

int main(int argc, char **argv)
{
	double tau,nu,v,lambdaf;
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE	 
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
		newfile >> tau >> lambdaf >> lambdat >> nu; // 8.6349e-9 0.1689e2 0.0360e2 2.637
    	} else
    	{
		exit(EXIT_FAILURE);
	}
    	newfile.close(); //close the file object.
    
    	//COMPUTING MU, ETA AND CRITICAL RADIUS
    	v=lambdaf/tau;
	mu = lambdat* v /3; 
	eta= v *( nu -1) / lambdaf;
	Numerics num;
	r1=num.bisection(CriticalRadius,0,10,1e-5);
	cout << endl;
	R=8.5;
	alpha=num.newton(NewAlpha,AlphaPrime, 0, 1e-8,100);
	
	//COMPUTING ns
	double n, dx,tfin,dt,A;
	double k = sqrt (( eta + alpha ) / mu ) ;
	cout << "This is now a small test for x and t; please write dx, final t and dt" << endl;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> tfin;
	cout << "dt: ";
	cin >> dt;
	int NPoints=2*abs(R/dx);
	int NStep=abs(tfin/dt);
	A=1/ (sin ( k * R ) / R);
	newfile.open("output.txt",ios::out);
	newfile << "n(r,t)" << "	" << "t" << "	" << "r" << endl;
	for (int c=1; c<NStep; ++c)
	{
		for (int i=1; i<NPoints; ++i)
		{
			n=A*exp(- alpha *c*dt)*sin(k*(-R+i*dx))/(-R+i*dx);
			newfile << n << "	" << c*dt << "	" << -R+i*dx << endl;	
		};
	};
	newfile.close();	
	return 0;
}
