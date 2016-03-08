#include "../include/smartuq.h"

using namespace smart;

int main(){

    /**
      17 uncertain variables (7 states + 10 parameters)
      4 degree polynomial expansion
      **/

    //variable allocation
    std::vector<chebyshev_polynomial<double> > x0, param;
    std::vector<chebyshev_polynomial<double> > xf;

    //initialise 7 state variables as Chebycheff base of order 1 in the variable i
    for(int i=0;i<7;i++){
        x0.push_back(chebyshev_polynomial<double>(17,4,i));
        x0[i].to_monomial_basis();
    }
    //initialise 10 parameters variables as Chebycheff base of order 1 in the variable 7+i
    for(int i=0;i<10;i++){
        param.push_back(chebyshev_polynomial<double>(17,4,i+7));
        param[i].to_monomial_basis();
    }

    x0[0].initialize_M(17,4);

    dynamics::twobody<chebyshev_polynomial<double> > dyn(param);
    integrator::rk4<chebyshev_polynomial<double> > integrator(&dyn);

    integrator.integrate(0,6000,100,x0,xf);

    x0[0].delete_M();
}
