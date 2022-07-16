#ifndef NUMERICS
#define NUMERICS

/* A simple class of useful numerical methods developed by Xotta Francesco
 * 
 */


 class Numerics
{
	public:
		double derivative (double (*func)(double), double dx, double pointx);
		double diffEq (double (*func) (double,double), double a, double b, double iniCond, double h);
		double secDifferentialEq (double (*func)(double, double), double (*func2)(double,double), double a, double b, double iniCond, double iniCondDer, double h);
		double integration(double (*func)(double), double a, double b, double h);
		double integrationSimpsRule(double (*func)(double), double a, double b, double h);
		double bisection(double (*f)(double), double a, double b, double error);
		double trapezoidal(double (*f)(double), double a, double b, double h);
		double trapezoidal(double (*f)(double, double), double ax, double bx, double ay, double by, double h);
		double trapezoidal(double (*f)(double, double, double), double ax, double bx, double ay, double by, double az, double bz, double h);
		double newton(double (*func)(double), double (*funcprime)(double), double p0, double err);
	private:
		double RK4Alg (double (*func) (double,double), double t, double y, double h);
};


#endif
