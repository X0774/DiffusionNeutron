#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"

using namespace std;

//DEFINITION OF VARIABLES
double L,A,mu,eta;
int N;
int p=0; //p,q,r must be on scope and must thus be defined before the code
int q=0;
int r=0;

double Function (double x, double y, double z)
{
	return (1 -( x / L ) ) *(1 -( y / L ) ) * (1 -( z / L ) ) * x * y * z /pow(( L /2),3); //The initial condition were different this time
};

double Int (double x, double y, double z) //It's way easier to define the integrand like this
{
	return 8*Function(x,y,z)*sin(p*M_PI*x/L)*sin(q*M_PI*y/L)*sin(r*M_PI*z/L)/(L*L*L);
}


int main(int argc, char **argv)
{
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
	  newfile >> L >> A >> N >> mu >> eta; //0.192 1 5 2.3446e+05 1.8958e+8
    	} else
    	{
		exit(EXIT_FAILURE);
	}
     	newfile.close(); //close the file object.
	
	
	////COMPUTE a[p][q][r]
	double*** a=new double** [N];
	for (int i=0;i<N;++i)
		{
			a[i]=new double* [N];
			for (int c=0; c<N; ++c)
				a[i][c]=new double [N];
		}		
	Numerics num;
	double n;
	for (p=0;p<N;++p)
	{
		for (q=0; q<N; ++q)
		{	
			for (r=0; r<N; ++r)
			{
				a[p][q][r]=num.trapezoidal(Int,0,L,0,L,0,L,0.01); 
				cout << a[p][q][r] << "	";
			};
			cout << endl;
		}
	};

	//WRITE n IN OUTPUT FILE
	cout << "This is now a small test for x and t; please write dx, computing time" << endl;
	double dx; //0.001
	double t; //2e-7;
	cout << "dx: ";
	cin >> dx;
	cout << "t: ";
	cin >> t;
	int NPoints=abs(L/dx); //it's a cubic box L*L*L
	newfile.open("output.txt",ios::out);
	newfile << "n(" << t <<",x,y,z)" << "	" << "x" << "	" << "y" << "	" << "z" <<endl;
	for (int i=1; i<NPoints; ++i)
		for (int k=1; k<NPoints; ++k)
			for (int l=1; l<NPoints; ++l)
			{
				n=0; //At the end of every loop in p,q,r, n must be reset
				for (p=0;p<N;++p)
				for (q=0; q<N;++q)
				for (r=0; r<N; ++r) //Indentation convention changed here for readability
					n=n+a[p][q][r] * exp( eta*(t)-mu*(  pow((p*M_PI/L),2)+pow((q*M_PI/L),2)+pow((r*M_PI/L),2) )*(t))  *  sin(p*M_PI*i*dx/L) * sin(q*M_PI*k*dx/L) * sin(r*M_PI*l*dx/L);
				newfile << n << "	" << 0+i*dx << "	" << 0+k*dx << "	" << 0+l*dx << endl;
			};	
	newfile.close();
	
	//DELETING DYNAMICALLY ALLOCATED POINTERS
	for (int i=0; i<N; ++i) 
		{
			for (int c=0; c<N; ++c)
				delete[] a[i][c];
			delete[] a[i];
		}
	delete[] a;
	return 0;
}
//Programmers swear to C++ pacific from the code 
