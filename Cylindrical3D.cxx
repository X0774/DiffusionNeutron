#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;

//DEFINITION OF VARIABLES NEEDED IN ALL SCOPES
double L,r1,A,mu,eta; 
int N;
double* alpha; //alpha[q] is the q-th zero of the bessel function of zeroth kind; we don't need to compute them each time we need them, so we are going to read them from an input file
int p=1; 
int q=0;

double Function (double r, double z)
{
	return  A*(1-pow(r/r1,2))*sin(M_PI*z/L);
};

double Int (double r, double z) //It's way easier to define the integrand like this
{
	return 4*   Function(r,z)*r* sin(M_PI*z*p/L) * j0(r*alpha[q]/r1)  /(L*pow(j1(alpha[q])*r1,2));
}


int main(int argc, char **argv)
{
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE	 
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
		newfile >> L >> r1 >> A >> N >> mu >> eta; //0.192 0.104 1 11 2.3446e+05 1.8958e+8
   	} else
    	{
		exit(EXIT_FAILURE);
	}
    	newfile.close(); //close the file object.

	alpha=new double[N];
	newfile.open("besselzeros.txt",ios::in);
	for (int i=0; i<N; ++i)
		{
			newfile >> alpha[i];
		};
	newfile.close();
	
	////COMPUTE a[p][q][r]
	double* a=new double [N]; //We don't need the dependency on p, we can show this later	
	Numerics num;
	double n;
	for (q=0;q<N;++q)
	{
		a[q]=num.trapezoidal(Int,0,r1,0,L,0.001); 
		cout << a[q] << "	";
		cout << endl;
	};

	//WRITE n IN OUTPUT FILE
	cout << "This is now a small test for x and t; please write dx, computing time" << endl;
	double dx; //0.001
	double t; //2e-7;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> t;
	int NPoints=abs(L/dx); 
	int NPointsr=abs(r1/dx);
	newfile.open("output.txt",ios::out);
	newfile << "n(" << t <<",r,z)" << "	" << "r" << "	" << "z" <<endl;
	for (int i=1; i<2*NPointsr; ++i)
		for (int k=1; k<NPoints; ++k)
		{
			n=0; //At the end of every loop in p,q,r, n must be reset
			for (q=0;q<N;++q)
				n=n+a[q] * j0((-r1+i*dx)*alpha[q]/r1) * sin(p*M_PI*k*dx/L) * exp( t * ((eta*pow(r1*L,2)-mu*( pow(alpha[q]*L,2) + pow(M_PI*p*r1,2)  )) )/pow(r1*L,2));
			newfile << n << "	" << -r1+i*dx << "	" << 0+k*dx << endl;
			};	
	newfile.close();
	
	//DELETING DYNAMICALLY ALLOCATED POINTERS
	delete[] a;
	delete[] alpha;
	return 0;
}


