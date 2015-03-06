#include <iostream>
#include <vector>

#include "chebyshev_polynomial.h"


#include <iterator>
#include <algorithm>


using namespace std;


int main()
{
    cout << "Welcome to Chebyshev Algebra!" << endl;

    Chebyshev_Polynomial<double> poly(3,5);

    std::cout<<poly<<std::endl;

    return 0;
}
