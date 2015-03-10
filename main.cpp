#include <iostream>
#include <vector>

#include "chebyshev_polynomial.h"


#include <iterator>
#include <algorithm>


using namespace std;


int main()
{
    cout << "Welcome to Chebyshev Algebra!" << endl;

    Chebyshev_Polynomial<double> poly1(2,3), poly2(2,3);
    double coeff1_data[] = {1,2,-1.0/3.0,0,3,0,0,0,0,0};
    double coeff2_data[] = {2,0,0,1,0,0,-5.0,0,1.0/2.0,0};
    std::vector<double> coeff1(coeff1_data, coeff1_data + sizeof(coeff1_data) / sizeof(double)), coeff2(coeff2_data, coeff2_data + sizeof(coeff2_data) / sizeof(double));

    poly1.set_coeffs(coeff1);
    poly2.set_coeffs(coeff2);

    std::cout<<poly1<<std::endl;
    std::cout<<poly2<<std::endl;

    Chebyshev_Polynomial<double> res = poly1*poly2;

    std::cout<<res<<std::endl;
    return 0;
}
