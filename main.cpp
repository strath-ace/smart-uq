#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "chebyshev_polynomial.h"
#include "elementary_functions.h"
#include "integrators.h"
#include "f.h"


#include <iterator>
#include <algorithm>


using namespace std;

int main()
{
    cout << "Welcome to Chebyshev Algebra!" << endl;

//    Chebyshev_Polynomial<double> poly1(2,3), poly2(2,3);
//    double coeff1_data[] = {1,2,-1.0/3.0,0,3,0,0,0,0,0};
//    double coeff2_data[] = {2,0,0,1,0,0,-5.0,0,1.0/2.0,0};
//    std::vector<double> coeff1(coeff1_data, coeff1_data + sizeof(coeff1_data) / sizeof(double)), coeff2(coeff2_data, coeff2_data + sizeof(coeff2_data) / sizeof(double));

//    poly1.set_coeffs(coeff1);
//    poly2.set_coeffs(coeff2);

//    std::cout<<poly1<<std::endl;
//    std::cout<<poly2<<std::endl;

//    Chebyshev_Polynomial<double> res = poly1*poly2;

      //Chebyshev_Polynomial<double> x(1,50,0);

      //NB: use double constants!!!
      //Chebyshev_Polynomial<double> f = x+7.0;//(4.0-x)*(4.0-x)*(5.0+x);
      //Chebyshev_Polynomial<double> f = x;

      //std::cout<<"f"<<std::endl;
      //std::cout<<f<<std::endl;

      //Chebyshev_Polynomial<double> g1 = 1.0/f;
      //Chebyshev_Polynomial<double> g2 = inv(f);
      //Chebyshev_Polynomial<double> g2 = sin(f);


      //std::cout<<g1<<std::endl;
      //std::cout<<g2<<std::endl;

    std::ofstream file;
    file.open ("results.out");

    std::vector<Chebyshev_Polynomial<double> > x0;
    x0.push_back(Chebyshev_Polynomial<double>(2,10));
    x0.push_back(Chebyshev_Polynomial<double>(2,10));
//    x0.push_back(Chebyshev_Polynomial<double>(4,10));
//    x0.push_back(Chebyshev_Polynomial<double>(4,10));

    std::cout<<x0[0]<<std::endl;

    x0[0].set_coeffs(1,1);
    x0[1].set_coeffs(2,1);
//    x0[2].set_coeffs(3,1);
//    x0[3].set_coeffs(4,1);

    std::vector<Chebyshev_Polynomial<double> > res;
    res.push_back(Chebyshev_Polynomial<double>(2,10));
    res.push_back(Chebyshev_Polynomial<double>(2,10));
//    res.push_back(Chebyshev_Polynomial<double>(4,10));
//    res.push_back(Chebyshev_Polynomial<double>(4,10));

//    //translation [a,b] --> [-1,1]
//    double a = 10.0;
//    double b = 13.0;
//    x0[0] = (2.0*x0[0] - (a+b))/(b-a);

    res = rk4<double>(f,x0,0.001,500);

    std::cout<<res[0];
    for(int i=0; i<res[0].get_coeffs().size(); i++){
        file<<left<<setw(16)<<res[0].get_coeffs()[i] <<left<<setw(16)<< res[1].get_coeffs()[i]<<"\n";
        //<<left<<setw(16)<<res[2].get_coeffs()[i] <<left<<setw(16)<< res[3].get_coeffs()[i]<<"\n";
    }
    file.close();


    return 0;
}
