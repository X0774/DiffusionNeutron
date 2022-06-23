#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include <string>
#include "Numerics.h"

using namespace std;
//DEFINITION OF VARIABLES
double L;
double A;
double lambda;
int N;
double mu;
double eta;
double alpha=10; //This is a test value

int p=0; //since p is used outside of the main we need to predefine p in order for it to be in every function's scope

double T (double t, double x)
{
	return x*(eta-alpha);
};

double X (double t, double x)
{
	return x*alpha;
};

double Xder(double t, double x)
{
	return (eta-alpha);
};

double Function (double t, double x)
{
	return A*exp( -lambda*(((x-L/2)/(L/2)))*((x-L/2)/(L/2)) );
};

double simpleShooting (double (*func)(double,double), double (*funcDer)(double,double), double a, double b, double conta, double contb, double dx, double alpha, double err){//Second derivatives solving, alpha is first derivative in a
	//y''=func->(y')'=func con condizione iniziale data dal guess che è alpha
	Numerics c;
	double beta=c.secDifferentialEq(func,a,b,conta,alpha,dx); //This gives me the function in b by solving a differential equation of second order with initial conditions y(a), alpha
	double beta1=c.secDifferentialEq(func,a,b,0,1,dx); //IMPORTANT: for the problem in question (SCHROEDINGER EQUATION) this works, but that's because of how the equation is made. DO NOT USE THIS METHOD ELSEWHERE
	//cout << beta<<','<<beta1<<endl; //REMOVE COMMENT TO TEST
	
	if (abs(beta-contb)>err){ 
		double newalpha=alpha-(beta-contb)/beta1;
		return simpleShooting(func,funcDer,a,b,conta,contb,dx,newalpha,err);
	}else
		return alpha;
}

int main(int argc, char **argv)
{
	
	
	//Reading Parameters from a Configuration File; the configuration is given as L A lambda N mu eta as will (hopefully) be written in the README
	fstream newfile;
	newfile.open("config.txt",ios::in);
	if (newfile.is_open())
	{
	  newfile >> L >> A >> lambda >> N >> mu >> eta;
      
    } else
    {
		exit(EXIT_FAILURE);
	}
      newfile.close(); //close the file object.
	
	Numerics s;
	double n;
	alpha=simpleShooting(X,Xder,0,L,0,0,0.01,alpha,0.1);
	cout << alpha;
	
	double dx=0.001;
	int NPoints=abs(L/dx);
	double tfin=0.00002;
	double dt=0.0000001;
	int NStep=abs(tfin/dt);
	newfile.open("output.txt",ios::out);
	newfile << "n(x,t)" << "	" << "t" << "	" << "x" << endl;
	for (int c=1; c<NStep; ++c)
	{
		for (int i=1; i<NPoints; ++i)
		{
			n=s.diffEq(T,0,c*dt,Function(0,i*dx),0.01)*s.secDifferentialEq(X,0,i*dx,0,eta-alpha,0.01);
			cout << s.diffEq(T,0,c*dt,Function(0,i*dx),0.01) << "	" << Function(0,i*dx) << endl;
			newfile << n << "	" << c*dt << "	" << 0+i*dx << endl;
		}
	};
	newfile.close();
	
	
	double l;
	cin >> l;
	return 0;
}
