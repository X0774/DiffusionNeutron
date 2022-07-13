#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;
//DEFINITION OF VARIABLES
double L,A,lambda,mu,eta;
int N;
int p=0; //since p,q are used outside of the main we need to predefine p,q in order for them to be in every function's scope
int q=0;

double Function (double x, double y)
{
	return (1 -( x / L ) ) *(1 -( y / L ) ) * x * y /pow(( L /4),2); //The initial condition were different this time
};

double Int (double x, double y) //It's way easier to define the integrand like this
{
	return 4*Function(x,y)*sin(p*M_PI*x/L)*sin(q*M_PI*y/L)/(L*L);
}


int main(int argc, char **argv)
{
	//Reading Parameters from a Configuration File
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
	  newfile >> L >> A >> lambda >> N >> mu >> eta; //15.7 1 100 6 2.3446e+05 1.8958e+8
    	} else
    	{
		exit(EXIT_FAILURE);
	}
      	newfile.close(); //close the file object.
	
	//GENERATING a[p][q]
	double** a=new double* [N];
	for (int i=0;i<N;++i)
		a[i]=new double[N];		
			
	Numerics num;
	double n;
	for (p=0;p<N;++p)
	{
		for (q=0; q<N; ++q)
		{
			a[p][q]=num.trapezoidal(Int,0,L,0,L,0.01); 
			cout << a[p][q] << "	";
		}
		cout << endl;
	}

	//WRITE NS IN OUTPUT FILE
	cout << "This is now a small test for x and t; please write dx, computing time" << endl;
	double dx; //0.001
	double t; //0.00001;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> t;
	int NPoints=abs(L/dx); //it's a square L*L
	newfile.open("output.txt",ios::out);
	newfile << "n(" << t <<",x,y)" << "	" << "x" << "	" << "y" << endl;
	for (int i=1; i<NPoints; ++i)
		for (int k=1; k<NPoints; ++k)
		{
			n=0; //At the end of every loop in p,q, n must be reset
			for (p=1;p<N;++p)
				for (q=1; q<N;++q)
				{
					n=n+a[p][q] * exp( eta*(t)-mu*(  pow((p*M_PI/L),2)+pow((q*M_PI/L),2)  )*(t))  *  sin(p*M_PI*i*dx/L) * sin(q*M_PI*k*dx/L);
				}
			newfile << n << "	" << 0+i*dx << "	" << 0+k*dx << endl;
		};	
	newfile.close();
	
	//DELETING DYNAMICALLY ALLOCATED POINTERS
	for (int i=0; i<N; ++i) 
		delete[] a[i];
	delete[] a;
	return 0;
}
