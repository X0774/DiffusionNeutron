#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"
using namespace std;

double eta,mu,lambdat,r1,R;
double A=1; //it serves no purpose but existing
double alpha=0; //If we want criticality
double k=0;

double Rfunc (double r){
	return (A*sin(k*r))/r; 
};

double Rfuncprime (double r){
	Numerics num;
	return num.derivative(Rfunc,0.00001,r);
};

double neumannBC(double r)
{
	return Rfuncprime(r)+3*Rfunc(r)/(2*lambdat);
};

int main(int argc, char **argv)
{
	double tau,nu,v,lambdaf;
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE	 
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
		newfile >> tau >> lambdaf >> lambdat >> nu; // 8.6349e-9 0.1689*100 0.0360*100 2.637
    	} else
    	{
		exit(EXIT_FAILURE);
	}
    	newfile.close(); //close the file object.
    
   	//COMPUTING MU, ETA AND CRITICAL RADIUS
    	v=lambdaf/tau;
	mu = lambdat* v /3; 
	eta= v *( nu -1) / lambdaf;
	k=sqrt((eta+alpha)/mu);
	Numerics num;
	//r1=num.bisection(criticalradius,0,10,1e-5);
	r1=num.bisection(neumannBC,0,10,1e-5);
	cout <<r1<< endl;
	return 0;
}
