#include <iostream>
#include "Numerics.h"
#include <cmath> 
#include <math.h> 

using namespace std;



double Numerics::derivative (double (*func)(double), double dx, double pointx){
	/**A simple 5 points formula that gives the derivative of func in pointx using dx as small distance*/
	return (func(pointx-2*dx)-8*func(pointx-dx)+8*func(pointx+dx)-func(pointx+2*dx))/(12*dx);
};

double Numerics::RK4Alg (double (*func) (double,double), double t, double y, double h)
{
	/**Simple algorithm made for the RK4 differential equation*/
	double k1,k2,k3,k4;
	k1=func(t,y);
	k2=func(t+h/2,y+h*k1/2);
	k3=func(t+h/2,y+h*k2/2);
	k4=func(t+h,y+h*k3);
	return y+h*(k1+2*k2+2*k3+k4)/6;
};

double Numerics::diffEq (double (*func) (double,double), double a, double b, double iniCond, double h) 
{
	/**RK4 Differential Equation Solver; dy/dt=func(t,y) from a to b with y(a)=iniCond with h=distance in each step*/
	int Nstep=abs(b-a)/h;
	
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
	/**Simple second order differential equation solver, with y''=func2*y'+func(t,y) from a to b, with y(a)=iniCond, y'(a)=iniCondDer with h=distance in each step*/
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
	/**Simpson 3/8 rule for integration; fast, but not very precise. It integrates func from a to b with a step distance of h*/
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

double Numerics::integrationSimpsRule(double (*func)(double), double a, double b, double h) //TO ADD double t0, (*func)(double, double)
{
	/**Simpson 1/3 rule for integration. It integrates func from a to b with a step distance of h*/
	int Nstep=abs(b-a)/(2*h);
	
	if (b<a)
	{
		double z=b;
		b=a;
		a=z;
	}

	double sum=func(a)+func(b);
	
	for (int i=1;i<Nstep;++i)
	{
		sum=sum+2*func(a+2*i*h)+4*func(a+(2*i-1)*h);
	};
	sum=sum-2*func(a+2*Nstep*h);
	
	return h*sum*1/3;
}

double Numerics::bisection(double (*f)(double), double a, double b, double error)
{
	/**Simple Root Finding algorythm that uses bisection; it solves between a and b; if f has more than one 0 between a and b it may not give any result*/
	double newpar=a+(b-a)/2;
	
	if (b<a)
	{
		double z=b;
		b=a;
		a=z;
	}
	
	if (abs(b-a)<abs(error) || abs(f(newpar))==0)
		return newpar;
		
	if (f(a)*f(newpar)<=0) //Recursive condition
		return bisection(f,a,newpar,error);
	else if (f(newpar)*f(b)<=0)
		return bisection(f,newpar,b,error);
	else
		{
			cout << "Error: number of zeros in the selection is not 1"<<endl;
			return a;
		}	
};

double Numerics::newton(double (*func)(double), double (*funcprime)(double), double p0, double err)
{
	/**Simple Root Finding algorythm that uses Newton method; it requires a function func and its derivative funcprime. It also wants a starting point p0*/
	double p=p0-func(p0)/funcprime(p0);
	if ((func(p) == 0) || (abs(p-p0) < err))
		return p;
	else 
		return newton(func,funcprime,p,err);
};

double Numerics::trapezoidal(double (*f)(double), double a, double b, double h)
{	
	/**Integral solver that uses the trapezoidal method; it solves between a and b with step distance h*/
	int Nstep=abs((b-a)/h);
	double Area=0;
	for (int i=0;i<Nstep;++i){
		Area=Area+(f(a+(i+1)*h)-f(a+i*h))*h/2+f(a+i*h)*h;
	}
	return Area;
};

double Numerics::trapezoidal(double (*f)(double, double), double ax, double bx, double ay, double by, double h)
{
	/**Integral solver that uses the trapezoidal method; it solves between [ax,bx] and [ay,by] with step distance h; it solves for y for each step of x*/
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
	/**Integral solver that uses the trapezoidal method; it solves between [ax,bx], [ay,by], [az,bz] with step distance h; it solves for z for each step of y, which solves for each step of x*/
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
