/*
 *i want to preface this code by saying that the original Maple code does not give the same result as this code, there is no way it can
 * 1) The a[p] are too small of a factor of 1/2; I solved the integral with various tools and observed this
 * Actually the paper contains a simplified formula for a[p] which agrees with my results :\
 */



#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;

//DEFINITION OF VARIABLES
double r1,mu,eta,A,lambda; 
int N;
int p=0; 

double Function (double r)
{
	return  1-pow(r/r1,2); // 
	//return A * exp ( - lambda * pow(( r / r1 ),2)) ; //if we want gaussian initial conditions
};

double Int (double r) //It's way easier to define the integrand like this
{
	return 2*Function(r) *r* sin(p*M_PI*r/r1) / r1; //we get 2 times the results given in the paper, which seems to be wrong
	//(2/ r1 ) * r * f ( r ) * sin (( p ) * Pi * r / r1  
}


int main(int argc, char **argv)
{
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE	 
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
		newfile >> r1 >> A>> lambda >> N >> mu >> eta; //0.115 1 100 31 2.3446e+05 1.8958e+8; A and lambda are if you want to use Gaussian
    	} else
    	{
		exit(EXIT_FAILURE);
	}
    	newfile.close(); //close the file object.
	
	////COMPUTE a[p][q][r]
	double* a=new double [N]; 
	Numerics num;
	for (p=0;p<N;++p)
	{
		a[p]=num.trapezoidal(Int,0,r1,0.01); 
		cout << a[p] << "	";
		cout << endl;
	};

	//WRITE n IN OUTPUT FILE
	double n;
	cout << "This is now a small test for x and t; please write dx, final time t and time step dt" << endl;
	double dx; //0.001
	double t, dt; //2e-7
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> t;
	cout << "dt: ";
	cin >> dt;
	int Nstep=abs(t/dt); 
	int NPointsr=abs(r1/dx);
	newfile.open("output.txt",ios::out);
	newfile << "n(t,r)" << "	" << "t" << "	" << "r" <<endl;
	for (int i=1; i<2*NPointsr; ++i)
	{	
		for (int k=1; k<Nstep; ++k)
		{
			 n=0;//At the end of every loop in p,q,r, n must be reset
			//if (-r1+i*dx!=0) If we want to avoid infinities in the denominator this is a rather expensive method
			{
				for (p=1;p<N;++p)
					//n=n+a[p] * exp((eta-mu*pow(p*M_PI/r1,2))*k*dt) *sin(p*M_PI*(-r1+i*dx)/r1)/(-r1+i*dx);
					n=n+ (a [ p ]/ (-r1+i*dx) ) * exp ( eta *k*dt - mu *pow((( p ) * M_PI / r1 ),2) *k* dt ) * sin (( p ) * M_PI * (-r1+i*dx) / r1 ); 
				newfile << n << "	" << 0+k*dt << "	" << -r1+i*dx << endl;
			};
		};	
	};
	newfile.close();
	
	//DELETING DYNAMICALLY ALLOCATED POINTERS
	delete[] a;
	return 0;
}
//Programmers swear to C++ pacific from the code 
