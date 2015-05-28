#ifndef F_H
#define F_H

#include <vector>
#include "chebyshev_polynomial.h"

//autonomous ODE system definition
template <class T>
std::vector<Chebyshev_Polynomial<T> > f(std::vector<Chebyshev_Polynomial<T> > x) {
	int n = x.size();
	int degree = x[0].get_degree();
	int nvar = x[0].get_nvar();
	std::vector<Chebyshev_Polynomial<T> > deriv;

	//allocate memory for polynomials
	for(int i=0; i<n; i++){
	    deriv.push_back(Chebyshev_Polynomial<T>(nvar,degree));
	}

// pendolo
	deriv[0] = x[1];
	deriv[1] = -9.81/2.0*sin(x[0]);

//2 body problems
//	deriv[0] = x[2];			   //r         --> x[0]
//	deriv[1] = x[3];			   //theta     --> x[1]
//	deriv[2] = x[0]*x[3]*x[3]-1.0/(x[0]*x[0]); //r_dot     --> x[2]
//	deriv[3] = -2.0*x[2]*x[3]/x[0];		   //theta_dot --> x[3]

//van der pol
//	deriv[0] = x[1];
//	deriv[1] = 0.5*(1.0-x[0]*x[0])*x[1]-x[0];

//logistic map
//	deriv[0] = x[1]*x[0]*(1.0-x[0]);
//	deriv[1] = x[1];

//lotka volterra
//        deriv[0] = 3.0*x[0]*(1.0-x[1]);
//        deriv[1] = 3.0*x[1]*(x[0]-1.0);

//chiara stuff :)
//        deriv[0] = x[1];
//        deriv[1] = -1.0*x[2]*x[0] - 1.0*x[1];
//        deriv[2] = 0.0;

//trivial
	  //deriv[0] = x[1];
	  //deriv[0] = x[0]*x[0];


	return deriv;
}


#endif // F_H
