# ReadME
Good morning, this is the first ReadME I have ever written, so please don't be mad if it's useless!

This project is set of codes meant to reproduce the results obtained in the paper "Neutron Diffusion" (https://www.researchgate.net/publication/323035158_Neutron_diffusion) by Graham W Griffiths; the codes roughly follow the Maple Examples given in said paper, solving the problem in the same way the author did.
I have personally created a Numerics class (keep in mind that, while this is not an optimal choice in a world where SymbolicC++ and similar lybraries exist, I wanted to develop those methods by myself to fully comprehend them).

# WHAT ARE THE FILES INSIDE THIS PROJECT:
There are 2 txt files (ignoring besselzeros, which was used to write the output of "dumbZeros.cxx", which was used to obtain Bessel's Functions zeros), which are config and output.
Config.txt contains all the initial datas that the programs must be given to work properly; suggested values (which follows from Table 7 of the aforementioned paper) are given as comments next to the lines in which config.txt is being read.
Output.txt serves the purpose of writing the output of the codes, which can then be read with GNUplot, Python or Matlab (C++ does not have a native way to plot functions); I personally used Matlab since UniBO allows its students to use it for free.

There are also 8 .cxx files (without considering dumbZeros.cxx), which are:
TEST1D.cxx: a small code meant to reproduce the neutron density of a 1-dimensional string of fissile material with L=L(critical); it produces an output of n(t,x)

TEST2D.cxx: same as TEST1D.cxx, but in 2 dimensions; it produces an output of n(t,x,y) with t fixed

Cartesian3D.cxx: same as TEST2D and TEST1D; it produces an output of n(t,x,y,z)

Cylindrical3D.cxx: code meant to reproduce the neutron density of a 3-dimensional cylinder of fissile material; it produces an output of n(t,r,z) with t fixed

Spherical3D.cxx: code meant to reproduce the neutron density of a 3-dimensional sphere of fissile material; it produces an output of n(t,r). it must be noted that I obtained results slightly different from the one obtained in the paper

SPHERICALNEUMAN.cxx: code meant to reproduce a sphere of fissile material with Neuman Boundary Conditions; it produces an output of n(t,r).
Notice that there is another file (alternativer1computation) which shows how one can compute r1 by solving Neumann Boundary Condition's equation (for alpha=1, so assuming criticality)

MEANFREEPATH.cxx: a small code meant to compute the mean free path of neutron escaping or reacting in the core; the output is the probability density function of neutron initiating a reaction in the block

The class Numerics.cxx is a class whose purpose is to compute as fast as possible integrals, derivatives and differential equations, working with standard C++ functions (keep in mind that this is not optimal, but was done just for clarity).

# INSTALLATION:
To use this code you need to compile it, something that can be done by using the command

  "g++ EXAMPLE.cxx Numerics.cxx"
 
Since the code is simple there isn't a dedicated directory for Numerics.cxx and Numerics.h; if you decided to create a separated directory you will need to:

1) Modify the code so that #include "Numerics.h" becomes #include "EXAMPLEDIRECTORY/Numerics.h"
2) Call from Command Prompt "g++ EXAMPLE.cxx EXAMPLEDIRECTORY/Numerics.cxx"

With the flag -Wall you can also experiment the absence of warnings!
# FAQ:
**Why did you use functions and not lambdas?**

Functions are more intuitive

**Why did you use arrays and not vectors**

I feel like arrays are more intuitive and are generally faster

**Why did you use trapezoidal instead of the integration method, which is faster?**

The integration method was developed using 3/8 rule and is thus a little bit less precise; I didn't need speed since the code worked pretty fast, so I kept using the Trapezoidal method

**Why are the bessel function zeros obtained with outside the code of Cylindrical3D.cxx?**

I thought it would have been way easier to compute Bessel Functions zeros externally and save them in a vector/array that could be read by the code instead of finding them from scratch each time I needed them

**What naming convention did you use?**

Methods have a lower case initial letter, functions have an upper case initial letter; dumbZero and SPHERICALNEUMAN are execptions (dumbZero because it's just a sample code that should not really be considered part of the project, SPHERICALNEUMAN because I didn't like having "Alpha" instead of "alpha"). MEANFREEPATH also has a function called prx, which comes from the name given in the original paper.

**What indentation convention did you use?**

I used tabs for every new scope
