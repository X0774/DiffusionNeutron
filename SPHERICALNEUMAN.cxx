/*
 * Symbolic computations are not supported and thus I will not follow some of the computations; I will also need to "cheat" the system in order to obtain some of the results
 * that could have been obtained with Mathematica, Matlab or Maple
 * 
 * For C++ we have SymbolicC++, which is a very good way to implement symbolic computations (but, since the numerical part is literally what the project is about, it would be cheating!)
 */ 

#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"
using namespace std;

double alphamin,alphamax,alpha,eta,mu,lambdat,r1,R;

double criticalradius(double r)
{
	return -1 + 3*r/lambdat/2 + r*sqrt(eta/mu)*cos(r*sqrt(eta/mu))/sin(r*sqrt(eta/mu));
};

double newalpha (double a)
{
	return -1 + 3*R/lambdat/2 + R*sqrt((eta+a)/mu)*cos(R*sqrt((eta+a)/mu))/sin(R*sqrt((eta+a)/mu));
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
	Numerics num;
	r1=num.Bisection(criticalradius,0,10,1e-3);
	R=8.5;
	alpha=num.Bisection(newalpha,-2,2e8,1e-3);
	cout << r1 <<"	"<< alpha;
	
	cin >> alpha;
	return 0;
}

