#include "main_list.h"
#include "utils.h"
void main_tests(){

    ////**********************************
    ////*******   TESTS EVALUATE   *******
    ////**********************************

    for (int nvar=4;nvar<=4;nvar++){
        cout << "[ ";
        for(int degree=4;degree<=12;degree++){
            Chebyshev_Polynomial<double> x (nvar, degree);
            cout << x.get_coeffs().size() << " , ";
        }
        cout <<" ]"<< endl;
    }
    // int nvar=4;
    // int degree=4;
    // Chebyshev_Polynomial<double> x (nvar, degree);

    // Chebyshev_Polynomial<double> p = (1+0.5*x+.3*x*x+4.6*x*x*x+18*x*x*x*x)/(x+2);

    // std::vector<double> point;
    // point.push_back(.6);

    // // cout << p << endl << endl;
    // cout << p.evaluate(point) << endl;
    
}