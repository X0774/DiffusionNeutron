#include <iostream>
#include <cmath> 
#include <math.h> 
#include <fstream>
#include "Numerics.h"
using namespace std;

//DEFINITION OF VARIABLES NEEDED IN ALL SCOPES
double R0, sigmaArea, s, n, sigma, NA, N0, rho, At;

double Prx (double x)
{
	return 1-exp(-sigma*n*x);
};

double prx (double x)
{
	Numerics num;
	return num.derivative(Prx,1e-3,x);
};

double Int (double x)
{
	return x*prx(x);
};

int main(int argc, char **argv)
{
	//Definition of variables
	double Prob1=0; //Probability that an individual neutron precipitates a reaction
	double Pesc=0; //Probability that a Neutron will escape
	double x=0;
	double Nesc=0; //Number of escaping neutrons
	double Nreact=0; //Number of reactions
	double z=0;
	double R=0; //Reaction Rate
	
	//READ INPUT FROM EXTERNAL CONFIGURATION FILE	 
	fstream newfile;
	newfile.open("config.txt",ios::in); 
	if (newfile.is_open())
	{
		newfile >> sigma >> NA>> rho >>At; //1.235e-28 6.022e23 18.71 235.04
    	} else
    	{
		exit(EXIT_FAILURE);
	}
    	newfile.close(); //close the file object.
    	n=( rho * NA/ At ) *pow(10,6) ;
	
	//Computing some quantities
	s=1e-5; //Must be small
	R0=1; //Whatever value is fine
	N0=1; //Whatever value is fine
	sigmaArea=1; //Whatever value is fine
	R=R0*sigmaArea*s*n*sigma;
	Prob1=R/(R0*sigmaArea); //Is this useless? Yes? I know.
	Pesc=1-Prob1;
	cout << "How thick will the slab be?"<<endl<<"x: ";
	cin>>x;
	cout << endl;
	Nesc=N0*pow(Pesc, x/s);
	z=s*n*sigma;
	cout << Nesc<<"	";
	Nesc=N0*pow(1+z,-sigma*n*x/z); //Nesc=N0*pow(pow(1+z,1/z),-sigma*n*x);
	cout << Nesc << "	"<<endl;
	Nesc=N0*exp(-sigma*n*x);
	Nreact=N0-Nesc;
	cout << "The number of reactions is: "<< Nreact << endl;
	
	//OUTPUT
	double dx=0;
	double xm=0;
	Numerics num;
	cout << "Please write x and dx, I will compute the probability density for the values of x given:" <<endl<<"x: ";
	cin>>x;
	cout<<"dx: ";
	cin>>dx;
	cout<<endl;
	newfile.open("output.txt",ios::out);
	newfile << "Prx(x)	" << "x" << endl;
	int NStep=abs(x/dx);
	xm=num.integration(Int,0,x,1e-3);
	cout<<"The expectation value is: "<<xm<<endl;
	for (int c=1; c<NStep; ++c)
	{
		newfile << prx(c*dx) <<"	"<< c*dx<<endl;	
	};
	newfile.close();
	return 0;
}
