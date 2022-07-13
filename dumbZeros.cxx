/******************************************************************************

THIS CODE IS THE DUMBEST POSSIBLE WAY THAT SOMEONE COULD USE TO FIND 
BESSEL FUNCTIONS' ZEROS, PLEASE REFRAIN FROM USING IT IN A SERIOUS ENVIRONMENT

*******************************************************************************/

#include <iostream>
#include <cmath>
#include "Numerics.h"
#include <fstream>
using namespace std;


void dumbZeros(double (*func)(double), double a, double dx, const int q, double error,double *array)
{
    	double zero=a;
    	int c;
    	Numerics num;
	for (int i=0; i<q; ++i)
	{
        c=0;
        do{
             c++;
        } while (func(zero+c*dx)*func(zero)>0);
        array[i]=num.bisection(func,zero,zero+c*dx,error);
        zero=zero+c*dx;
    };
};

int main()
{
    const int N=10;
    double array[N];
    dumbZeros(j0,0,0.00001,N,0.00001,array);
    fstream newfile;
    newfile.open("besselzeros.txt",ios::out);
    cout << "ciao";
    for (int i=0;i<N; ++i)
        newfile << array[i] << "	";
    newfile.close();
    return 0;
}

/*
░░░░░░▄▀▒▒▒▒░░░░▒▒▒▒▒▒▒▒▒▒▒▒▒█
░░░░░█▒▒▒▒░░░░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒█
░░░░█▒▒▄▀▀▀▀▀▄▄▒▒▒▒▒▒▒▒▒▄▄▀▀▀▀▀▀▄
░░▄▀▒▒▒▄█████▄▒█▒▒▒▒▒▒▒█▒▄█████▄▒█
░█▒▒▒▒▐██▄████▌▒█▒▒▒▒▒█▒▐██▄████▌▒█
▀▒▒▒▒▒▒▀█████▀▒▒█▒░▄▒▄█▒▒▀█████▀▒▒▒█
▒▒▐▒▒▒░░░░▒▒▒▒▒█▒░▒▒▀▒▒█▒▒▒▒▒▒▒▒▒▒▒▒█
▒▌▒▒▒░░░▒▒▒▒▒▄▀▒░▒▄█▄█▄▒▀▄▒▒▒▒▒▒▒▒▒▒▒▌
▒▌▒▒▒▒░▒▒▒▒▒▒▀▄▒▒█▌▌▌▌▌█▄▀▒▒▒▒▒▒▒▒▒▒▒▐
▒▐▒▒▒▒▒▒▒▒▒▒▒▒▒▌▒▒▀███▀▒▌▒▒▒▒▒▒▒▒▒▒▒▒▌
▀▀▄▒▒▒▒▒▒▒▒▒▒▒▌▒▒▒▒▒▒▒▒▒▐▒▒▒▒▒▒▒▒▒▒▒█
▀▄▒▀▄▒▒▒▒▒▒▒▒▐▒▒▒▒▒▒▒▒▒▄▄▄▄▒▒▒▒▒▒▄▄▀
▒▒▀▄▒▀▄▀▀▀▄▀▀▀▀▄▄▄▄▄▄▄▀░░░░▀▀▀▀▀▀
▒▒▒▒▀▄▐▒▒▒▒▒▒▒▒▒▒▒▒▒▐
 */
