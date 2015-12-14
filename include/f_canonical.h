#ifndef F_CANONICAL_H
#define F_CANONICAL_H

#include <vector>
#include "Polynomial/canonical.h"

//autonomous ODE system definition
template <class T>
std::vector<Canonical_Polynomial<T> > f(const std::vector<Canonical_Polynomial<T> > &x, const std::vector<Canonical_Polynomial<T> > &param) {
	int n = x.size();
	int degree = x[0].get_degree();
	int nvar = x[0].get_nvar();
	std::vector<Canonical_Polynomial<T> > deriv;

	//allocate memory for polynomials
	for(int i=0; i<n; i++){
	    deriv.push_back(Canonical_Polynomial<T>(nvar,degree));
	}

// pendolo
//	deriv[0] = x[1];
//	deriv[1] = -9.81/2.0*sin(x[0]);

//2 body problems
//	deriv[0] = x[2];			   //r         --> x[0]
//	deriv[1] = x[3];			   //theta     --> x[1]
//	deriv[2] = x[0]*x[3]*x[3]-1.0/(x[0]*x[0]); //r_dot     --> x[2]
//	deriv[3] = -2.0*x[2]*x[3]/x[0];		   //theta_dot --> x[3]

	//2 body problems
//	deriv[0] = x[2];			   //r         --> x[0]
//	deriv[1] = x[3];			   //theta     --> x[1]
//	deriv[2] = x[0]*x[3]*x[3]-1.0/(x[0]*x[0]) + x[4]*sin(x[1]); //r_dot     --> x[2]
//	deriv[3] = (-2.0*x[2]*x[3] + x[4]*cos(x[1]))/x[0];  //theta_dot --> x[3]
//	deriv[4] = 0;


// // //Accelerated Kepler problem planar
//  T mu = 1.0;
// 	Canonical_Polynomial<T> tmp_2D =  mu/pow(sqrt(x[0]*x[0]+x[1]*x[1]), 3);
// 	deriv[0] =  x[2];		// vx
// 	deriv[1] =  x[3];		// vy
// 	deriv[2] = -1.0*tmp_2D*x[0];//+ param[1];	// ax
//     deriv[3] = -1.0*tmp_2D*x[1] + param[0];	// ay

// //Accelerated Kepler problem spatial
	// T mu = 1.0;
	// Canonical_Polynomial<T> tmp_3D =  mu/pow(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), 3);
	// deriv[0] = x[3];
	// deriv[1] = x[4];
	// deriv[2] = x[5];
	// deriv[3] = -1.0*tmp_3D*x[0]+param[1];
	// deriv[4] = -1.0*tmp_3D*x[1]+param[0];
	// deriv[5] = -1.0*tmp_3D*x[2]+param[2];

// //van der pol
	// T mu = .5;
	// deriv[0] = x[1];
	// deriv[1] = mu*(1.0-x[0]*x[0])*x[1]-x[0];

// Lotka-Volterra
	Canonical_Polynomial<T> a=param[0];
	Canonical_Polynomial<T> b=param[1];
	Canonical_Polynomial<T> c=param[2];
	Canonical_Polynomial<T> d=param[3];
	
	Canonical_Polynomial<T> xy = x[0]*x[1];
	deriv[0] = a*x[0]-b*xy;
	deriv[1] = -c*x[1]+d*xy;

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
//	  deriv[0] = 1.0/x[0];


	return deriv;
}


#endif // F_CANONICAL_H
