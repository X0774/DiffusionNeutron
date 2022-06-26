Member Function Documentation
◆ bisection()
double Numerics::bisection	(	double(*)(double) 	f,
double 	a,
double 	b,
double 	error 
)		
Simple Root Finding algorythm that uses bisection; it solves between a and b; if f has more than one 0 between a and b it may not give any result

◆ derivative()
double Numerics::derivative	(	double(*)(double) 	func,
double 	dx,
double 	pointx 
)		
A simple 5 points formula that gives the derivative of func in pointx using dx as small distance

◆ diffEq()
double Numerics::diffEq	(	double(*)(double, double) 	func,
double 	a,
double 	b,
double 	iniCond,
double 	h 
)		
RK4 Differential Equation Solver; dy/dt=func(t,y) from a to b with y(a)=iniCond with h=distance in each step

◆ integration()
double Numerics::integration	(	double(*)(double) 	func,
double 	a,
double 	b,
double 	h 
)		
Simpson 3/8 rule for integration; fast, but not very precise. It integrates func from a to b with a step distance of h

◆ newton()
double Numerics::newton	(	double(*)(double) 	func,
double(*)(double) 	funcprime,
double 	p0,
double 	err 
)		
Simple Root Finding algorythm that uses Newton method; it requires a function func and its derivative funcprime. It also wants a starting point p0

◆ RK4Alg()
double Numerics::RK4Alg	(	double(*)(double, double) 	func,
double 	t,
double 	y,
double 	h 
)		
private
Simple algorithm made for the RK4 differential equation

◆ secDifferentialEq()
double Numerics::secDifferentialEq	(	double(*)(double, double) 	func,
double(*)(double, double) 	func2,
double 	a,
double 	b,
double 	iniCond,
double 	iniCondDer,
double 	h 
)		
Simple second order differential equation solver, with y''=func2*y'+func(t,y) from a to b, with y(a)=iniCond, y'(a)=iniCondDer with h=distance in each step

◆ trapezoidal() [1/3]
double Numerics::trapezoidal	(	double(*)(double) 	f,
double 	a,
double 	b,
double 	h 
)		
Integral solver that uses the trapezoidal method; it solves between a and b with step distance h

◆ trapezoidal() [2/3]
double Numerics::trapezoidal	(	double(*)(double, double) 	f,
double 	ax,
double 	bx,
double 	ay,
double 	by,
double 	h 
)		
Integral solver that uses the trapezoidal method; it solves between [ax,bx] and [ay,by] with step distance h; it solves for y for each step of x

◆ trapezoidal() [3/3]
double Numerics::trapezoidal	(	double(*)(double, double, double) 	f,
double 	ax,
double 	bx,
double 	ay,
double 	by,
double 	az,
double 	bz,
double 	h 
)		
Integral solver that uses the trapezoidal method; it solves between [ax,bx], [ay,by], [az,bz] with step distance h; it solves for z for each step of y, which solves for each step of x
