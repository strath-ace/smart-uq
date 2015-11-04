#include "main_list.h"
#include "utils.h"

void main_tests(){

    srand(time(NULL));

    // Canonical_Polynomial <double> x (1,5,0);
    // Canonical_Polynomial <double> p=1+x*x;
    // std::vector<Canonical_Polynomial<double> > q;
    // q.push_back(3+x);

    // cout << p << endl;
    // cout << q[0] << endl;

    // cout << p*p << endl;;

    // cout << p.composition(q) << endl;




    // int nvar=8;
    // int deg=4;

    // Chebyshev_Polynomial <double> p (nvar,deg,0);


    // for (int j=0; j<p.get_coeffs().size();j++){
    //     p.set_coeffs(j,(double) (rand()%1000)/50.0-10); //randomly
    // }

    // p+=1000000000;
    // cout << p << endl;
    // cout << p/p << endl;



    Canonical_Polynomial<double> x (3,3,0);
    Canonical_Polynomial<double> y (3,3,1);
    Canonical_Polynomial<double> z (3,3,2);
    Canonical_Polynomial<double> p (3,3);

    p = 3 + 2*x*x*x + x*y + y*y*z + 5*z*z;

    std::vector<double> point (3,0);
    point[0]=1;
    point[1]=-.5;

    cout << p.evaluate(point) << endl;


 
} 