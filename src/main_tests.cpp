#include "main_list.h"

void main_tests(){

    int nvar=3;
    int degree=100;
    Chebyshev_Polynomial<double> x (nvar, degree,0);
    Chebyshev_Polynomial<double> y (nvar, degree,1);
    Chebyshev_Polynomial<double> z (nvar, degree,2);

    Chebyshev_Polynomial<double> p = 3*x+y+z;

    std::vector<double> point;
    point.push_back(.6);
    point.push_back(-1);
    point.push_back(.44);

    cout << p << endl << endl;
    cout << p.evaluate(point) << endl;
}
