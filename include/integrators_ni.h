#ifndef INTEGRATORS_NI_H
#define INTEGRATORS_NI_H

#include <vector>
#include "Polynomial/chebyshev.h"

// // Euler forward integrator step
// // Inputs
// //   x         dependent variable (initial values)
// //   t         independent variable (initial time)
// //   h          fixed step size

template <class T>
std::vector<T> euler(std::vector<T> (*f)(const std::vector<T> &x, const std::vector<T> &param), std::vector<T> &x, const std::vector<T> &param, const T &h){

        int n = x.size();

        std::vector<T> k1(n,0.0);

        //* Evaluate k1 = f(x).
        k1 = (*f)(x,param);

        for(int i=0; i<n; i++)
            x[i] += h*k1[i];

        return x;
}

// template <class T>
// std::vector<Chebyshev_Polynomial<T> > discrete_method(std::vector<Chebyshev_Polynomial<T> > (*f)(const std::vector<Chebyshev_Polynomial<T> > &x, const std::vector<Chebyshev_Polynomial<T> > &param), std::vector<Chebyshev_Polynomial<T> > &x, const std::vector<Chebyshev_Polynomial<T> > &param, const T &h) {

//         int n = x.size();
//         int degree = x[0].get_degree();
//         int nvar = x[0].get_nvar();

//         std::vector<Chebyshev_Polynomial<T> > k1;
//         for(int i=0; i<n; i++){
//             k1.push_back(Chebyshev_Polynomial<T>(nvar,degree));
//         }

//         //* Evaluate k1 = f(x).
//         k1 = (*f)(x,param);

//         for(int i=0; i<n; i++)
//             x[i] = k1[i];

//         return x;
// }



// // Runge-Kutta integrator step (2th order)
// // Inputs
// //   x         dependent variable (initial values)
// //   t         independent variable (initial time)
// //   h          fixed step size
// //   param      parameters

// template <class T>
// std::vector<Chebyshev_Polynomial<T> > rk2(std::vector<Chebyshev_Polynomial<T> > (*f)(const std::vector<Chebyshev_Polynomial<T> > &x, const std::vector<Chebyshev_Polynomial<T> > &param), std::vector<Chebyshev_Polynomial<T> > &x, const std::vector<Chebyshev_Polynomial<T> > &param, const T &h) {

//         int n = x.size();
//         int degree = x[0].get_degree();
//         int nvar = x[0].get_nvar();

//         std::vector<Chebyshev_Polynomial<T> > k1;
//         std::vector<Chebyshev_Polynomial<T> > k2;
//         std::vector<Chebyshev_Polynomial<T> > xtemp;

//         for(int i=0; i<n; i++){
//             k1.push_back(Chebyshev_Polynomial<T>(nvar,degree));
//             k2.push_back(Chebyshev_Polynomial<T>(nvar,degree));
//             xtemp.push_back(Chebyshev_Polynomial<T>(nvar,degree));
//         }

//         //* Evaluate k1 = f(x).
//         k1 = (*f)(x,param);

// 	//* Evaluate k2 = f(x+h*k1),
// 	for(int i=0; i<n; i++)
// 	    xtemp[i] = x[i]+k1[i]*h;
// 	k2 = (*f)(xtemp,param);

//         //* Return x(t+h) computed from second-order Runge Kutta.
//         for(int i=0; i<n; i++)
//             x[i] += (k1[i]+k2[i])*h/2.0;

//         return x;

// }


// // Runge-Kutta integrator step (4th order)
// // Inputs
// //   x         dependent variable (initial values)
// //   t         independent variable (initial time)
// //   h          fixed step size

template <class T>
std::vector<T> rk4 (std::vector<T> (*f)(const std::vector<T> &x, const std::vector<T> &param), std::vector<T> &x, const std::vector<T> &param, const T &h) {
    int n = x.size();
    // int degree = x[0].get_degree();
    // int nvar = x[0].get_nvar();

    std::vector<T> k1(n,0.0);
    std::vector<T> k2(n,0.0);
    std::vector<T> k3(n,0.0);
    std::vector<T> k4(n,0.0);
    std::vector<T> xtemp(n,0.0);

    // for(int i=0; i<n; i++){
    //     k1.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    //     k2.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    //     k3.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    //     k4.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    //     xtemp.push_back(Chebyshev_Polynomial<T>(nvar,degree));
    // }

    //* Evaluate k1 = f(x).
    k1 = (*f)(x,param);

    //* Evaluate k2 = f(x+h/2*k1),
    for(int i=0; i<n; i++)
        xtemp[i] = x[i]+k1[i]*h/2.0;
    k2 = (*f)(xtemp,param);

    //* Evaluate k3 = f(x+h/2*k2),
    for(int i=0; i<n; i++)
        xtemp[i] = x[i]+k2[i]*h/2.0;
    k3 = (*f)(xtemp,param);

    //* Evaluate k4 = f(x+h*k3),
    for(int i=0; i<n; i++)
        xtemp[i] = x[i]+k3[i]*h;
    k4 = (*f)(xtemp,param);

    //* Return x(t+h) computed from second-order Runge Kutta.
    for(int i=0; i<n; i++)
        x[i] += (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*h/6.0;

    return x;
}

#endif // INTEGRATORS_NI_H
