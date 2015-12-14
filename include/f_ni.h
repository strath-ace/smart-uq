#ifndef F_NI_H
#define F_NI_H

#include <vector>

//autonomous ODE system definition
template <class T>
std::vector<T> f(const std::vector<T> &x, const std::vector<T> &param) {
	int n = x.size();
	// int degree = x[0].get_degree();
	// int nvar = x[0].get_nvar();
	std::vector<T> deriv (n,0.0);

	//allocate memory for polynomials
	// for(int i=0; i<n; i++){
	//     deriv.push_back(Chebyshev_Polynomial<T>(nvar,degree));
	// }

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


// //Accelerated Kepler problem planar
// 	T mu = 1.0;
// 	T tmp_2D =  mu/pow(sqrt(x[0]*x[0]+x[1]*x[1]), 3);
// 	deriv[0] =  x[2];		// vx
// 	deriv[1] =  x[3];		// vy
// 	deriv[2] = -1.0*tmp_2D*x[0];// + param[1];	// ax
//  	deriv[3] = -1.0*tmp_2D*x[1] + param[0];	// ay

// //Accelerated Kepler problem spatial
 	// T mu = 1.0;
	// T tmp_3D =  mu/pow(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), 3);
	// deriv[0] = x[3];
	// deriv[1] = x[4];
	// deriv[2] = x[5];
	// deriv[3] = -1.0*tmp_3D*x[0]+param[1];
	// deriv[4] = -1.0*tmp_3D*x[1]+param[0];
	// deriv[5] = -1.0*tmp_3D*x[2]+param[2];

// van der pol
	// T mu = 0.5;
	// deriv[0] = x[1];
	// deriv[1] = mu*(1.0-x[0]*x[0])*x[1]-x[0];


// Lotka-Volterra
	T a, b, c, d;
	a = 1.0;//param[0];
	b = 1.0;//param[1];
	c = 1.0;//param[2];
	d = 1.0;//param[3];
	
	T xy = x[0]*x[1];
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


#endif // F_NI_H
