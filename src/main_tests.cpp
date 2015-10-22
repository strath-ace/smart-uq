#include "main_list.h"

void main_tests(){

    int nvar=1;
    int degree=100;
    Chebyshev_Polynomial<double> x (nvar, degree,0);

    Chebyshev_Polynomial<double> p = (1+0.5*x+.3*x*x+4.6*x*x*x+18*x*x*x*x)/(x+2);

    std::vector<double> point;
    point.push_back(.6);

    // cout << p << endl << endl;
    cout << p.evaluate(point) << endl;
}
