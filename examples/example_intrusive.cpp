/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../include/smartuq.h"
#include <fstream>

using namespace smartuq;

int main(){

    /**
      17 uncertain variables (7 states + 10 parameters)
      4 degree polynomial expansion
      **/

    int nvar = 7;
    int nparam  = 10;
    int poly_degree = 4;

    //polynomial allocation
    std::vector<chebyshev_polynomial<double> > x0, param;
    std::vector<chebyshev_polynomial<double> > xf;

    //initialisation ranges and constants terms
    double sma = 7378*pow(10,3); //*********
    double period = 2.0*M_PI/pow(sma,-3.0/2.0)/sqrt(398600.4415*pow(10,9));
    double tstart = 0;
    double tf = 0;
    double deltat = 0;
    double tend = period;

    std::vector<double> x(nvar), p(nparam), unc_x(nvar), unc_p(nparam);
    std::vector<std::vector<double> > ranges_x,ranges_p;

    x[0] = 7338*pow(10,3);
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 2*M_PI*sma/period;
    x[5] = 0.0;
    x[6] = 2000;

    p[0] = 0;
    p[1] = 0.5;
    p[2] = 0;
    p[3] = 1/30000;

    p[4] = 5.245*pow(10,-15);
    p[5] = 181050;
    p[6] = 4.4;

    p[7] = 0;
    p[8] = 0;
    p[9] = 0;

    // for (int i= 0 ; i < 7; i++) unc_x[i] = 0.000001;
    unc_x[0] = 1000;
    unc_x[1] = 1000;
    unc_x[2] = 1000;
    unc_x[3] = 5;
    unc_x[4] = 5;
    unc_x[5] = 5;
    unc_x[6] = 1.0;

    // for (int i= 0 ; i < 10; i++) unc_p[i] = 0.000001;
    unc_p[0] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[1] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[2] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[3] = 0.05*p[3];

    unc_p[4] = 0.01*p[4];
    unc_p[5] = 0.01*p[5];
    unc_p[6] = 0.01*p[6];

    unc_p[7] = 0.0001;
    unc_p[8] = 0.0001;
    unc_p[9] = 0.0001;

    for(int i=0; i<nvar; i++){
        ranges_x.push_back(std::vector<double>(2));
        ranges_x[i][0] = x[i]-unc_x[i];
        ranges_x[i][1] = x[i]+unc_x[i];
    }
    for(int i=0; i<nparam; i++){
        ranges_p.push_back(std::vector<double>(2));
        ranges_p[i][0] = p[i]-unc_p[i];
        ranges_p[i][1] = p[i]+unc_p[i];
    }

    //initialise 7 state variables as Chebycheff base of order 1 in the variable i
    for(int i=0;i<7;i++){
        x0.push_back(chebyshev_polynomial<double>(nvar+nparam,poly_degree,i,ranges_x[i][0], ranges_x[i][1], true));
    }
    //initialise 10 parameters variables as Chebycheff base of order 1 in the variable 7+i
    for(int i=0;i<10;i++){
        param.push_back(chebyshev_polynomial<double>(nvar+nparam,poly_degree,i+7,ranges_p[i][0], ranges_p[i][1], true));
    }

    clock_t begin, end;
    begin=clock();
    x0[0].initialize_M(nvar+nparam,poly_degree);
    std::vector<std::vector<double> > coeffs_all;

    dynamics::twobody<chebyshev_polynomial<double> > dyn(param);
    integrator::rk4<chebyshev_polynomial<double> > integrator(&dyn);

    deltat = 1000;
    for(int i=0; i+1<tend/deltat; i++){
        tf += deltat;
        integrator.integrate(tstart,tf,100,x0,xf);
        x0 = xf;
        tstart=tf;
        for(int j=0; j<nvar; j++){
            std::vector<double> coeffs = xf[j].get_coeffs();
            coeffs_all.push_back(coeffs);
        }
    }

    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "time elapsed : " << time_akp << endl << endl;
    x0[0].delete_M();

    std::ofstream file;
    file.open ("twobodyproblem.txt");
    for(unsigned int k=0; k<coeffs_all.size(); k++){
        for(unsigned int kk=0; kk<coeffs_all[k].size(); kk++)
            file << setprecision(16) << coeffs_all[k][kk] << " ";
        file << "\n";
    }
    file << "\n\n\n\n";
    file.close();
}
