#include <iostream>
#include "Numerics.h"
#include <cmath> 
#include <math.h> 

using namespace std;



double Numerics::derivative (double (*func)(double), double dx, double pointx){
	return (func(pointx-2*dx)-8*func(pointx-dx)+8*func(pointx+dx)-func(pointx+2*dx))/(12*dx);
};

double Numerics::RK4Alg (double (*func) (double,double), double t, double y, double h)
{
    double k1,k2,k3,k4;
    k1=func(t,y);
    k2=func(t+h/2,y+h*k1/2);
    k3=func(t+h/2,y+h*k2/2);
    k4=func(t+h,y+h*k3);
    return y+h*(k1+2*k2+2*k3+k4)/6;
};

double Numerics::diffEq (double (*func) (double,double), double a, double b, double iniCond, double h) //RK4 Differential Equation Solver; dy/dt=func(t,y)
{
    int Nstep=abs(b-a)/h;
    
    /*if (b<a)
    {
        double z=b;
        b=a;
        a=z;
    }*/
    
    double yi=iniCond;
    double yii;
    double t=a;
    for (int i=0; i<Nstep;++i)
        {
            yii=RK4Alg(func,t,yi,h);
            t=t+h;
        }
    return yii;    
};

double Numerics::secDifferentialEq (double (*func)(double, double), double (*func2)(double,double), double a, double b, double iniCond, double iniCondDer, double h){ //Simple second order differential Equation, y''=func2*y'+func(t,y)
	int Nstep=abs(b-a)/h;
	
	/*if (b<a)
    {
        double z=b;
        b=a;
        a=z;
    }*/
	
	double z=iniCond;
	double x=z+iniCondDer*h; //yii=yi+y'*dx	
	double yiii,yii,yi; //yi is j-1, yii is j, yiii is j+1
	for (int i=0;i<Nstep;++i){
		yi=z;
		yii=x;
		yiii=2*yii-yi+h*h*func(a+i*h,yii)+func2(a+i*h,yii)*(yii-yi)*h;
		x=yiii;
		z=yii;
	}
	return yiii;
};

double Numerics::integration(double (*func)(double), double a, double b, double h) //TO ADD double t0, (*func)(double, double)
{
    int Nstep=abs(b-a)/h;
    int N3=Nstep/3;
    
    if (b<a)
    {
        double z=b;
        b=a;
        a=z;
    }

    double sum=0;
    
    for (int i=0;i<N3;++i)
    {
        sum=sum+func(a+3*i*h)+3*func(a+3*i*h+h)+3*func(a+3*i*h+2*h)+func(a+3*i*h+3*h);
    }
    return h*sum*3/8;
    
}

double Numerics::Bisection(double (*f)(double), double a, double b, double error)
{
	double newpar=a+(b-a)/2;
	
	if (b<a)
    {
        double z=b;
        b=a;
        a=z;
    }
    
	if (abs(b-a)<abs(error) || abs(f(newpar))<abs(error))
		return newpar;
		
	if (f(a)*f(newpar)<=0) //Recursive condition
		return Bisection(f,a,newpar,error);
	else if (f(newpar)*f(b)<=0)
		return Bisection(f,newpar,b,error);
	else
		{
			cout << "Error: number of zeros in the selection is not 1";
			return a;
		}
		
};

double Numerics::Newton(double (*func)(double), double (*funcprime)(double), double p0, double err)
{
    double p=p0-func(p0)/funcprime(p0);
    if ((func(p) == 0) || (abs(p-p0) < err))
        return p;
    else 
        return Newton(func,funcprime,p,err);
};

double Numerics::trapezoidal(double (*f)(double), double a, double b, double h)
{	
	int Nstep=abs((b-a)/h);
	double Area=0;
	for (int i=0;i<Nstep;++i){
		Area=Area+(f(a+(i+1)*h)-f(a+i*h))*h/2+f(a+i*h)*h;
	}
	return Area;
};

double Numerics::trapezoidal(double (*f)(double, double), double ax, double bx, double ay, double by, double h)
{
	//We want to integrate by strings, so for constant x we integrate over y
	
	int Nstepx=abs((bx-ax)/h);
	int Nstepy=abs((by-ay)/h);	
	double Area=0;
	
	
	for (int i=0;i<Nstepx;++i)
	    for (int c=0; c<Nstepy; ++c)
	        Area=Area+(f(ax+(i+1)*h,ax+(c+1)*h)-f(ax+i*h,ay+c*h))*h*h/2+f(ax+i*h,ay+c*h)*h*h;
	return Area;
};

double Numerics::trapezoidal(double (*f)(double, double, double), double ax, double bx, double ay, double by, double az, double bz, double h)
{
	//We want to integrate by strings, so for constant x, y, we integrate over z
	
	int Nstepx=abs((bx-ax)/h);
	int Nstepy=abs((by-ay)/h);	
	int Nstepz=abs((bz-az)/h);	
	double Area=0;
	
	
	for (int i=0;i<Nstepx;++i)
	    for (int c=0; c<Nstepy; ++c)
			for (int k=0; k<Nstepz; ++k)
				Area=Area+(f(ax+(i+1)*h,ay+(c+1)*h,az+(k+1)*h)-f(ax+i*h,ay+c*h,az+k*h))*h*h*h/2+f(ax+i*h,ay+c*h,az+k*h)*h*h*h;
	return Area;
};
