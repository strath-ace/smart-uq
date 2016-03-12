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

int main(){

    /**
      17 uncertain variables (7 states + 10 parameters)
      4 degree polynomial expansion
      **/

    int nvar = 7;
    int nparam  = 10;
    int poly_degree = 4;

    int nsamples = combination(nvar+nparam,poly_degree);

    //polynomial allocation
    std::vector<double> x0, param;
    std::vector<double> xf;

    //initialisation ranges and constants terms
    double sma = 7378*pow(10,3); //*********
    double period = 2.0*M_PI/pow(sma,-3.0/2.0)/sqrt(398600.4415*pow(10,9));
    double tstart = 0;
    double tf = 0;
    double deltat = 0;
    double tend = period;

    std::vector<double> x(nvar), p(nparam), unc_x(nvar), unc_p(nparam);
    std::vector<double> ranges_lb(nvar+nparam), ranges_ub(nvar+nparam);

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
    p[3] = 1/30000.0;

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
        ranges_lb[i] = x[i]-unc_x[i];
        ranges_ub[i] = x[i]+unc_x[i];
    }
    for(int i=0; i<nparam; i++){
        ranges_lb[nvar+i] = p[i]-unc_p[i];
        ranges_ub[nvar+i] = p[i]+unc_p[i];
    }

    chebyshev_polynomial<double> poly(nvar+nparam,poly_degree, ranges_lb, ranges_ub);

    sampling::lhs<double> lhs_gen(nvar+nparam,nsamples,ranges_lb, ranges_ub);

    //initialise dynamics and integrator
    std::vector<std::vector<double> > coeffs_all;

    clock_t begin,end;
    begin=clock();

    // construct LHS sampling
    std::vector<std::vector<double> > LHS, H;
    for(int i=0;i<nsamples;i++){
        std::vector<double> sample=lhs_gen();
        LHS.push_back(sample);
    }

    deltat = 1000;
    for(int i=0; i<tend/deltat; i++){
        std::vector<std::vector<double> > y;
        std::vector<std::vector<double> > res_coeffs;
        tf += deltat;

        // perform nsamples forward integrations
        for(int i=0;i<nsamples;i++){
            std::vector<double> LHS_p, LHS_x;
            // separate states from parameters in the LHS sampling
            for(int j=0;j<nvar+nparam; j++){
                if(j<7)
                    LHS_x.push_back(LHS[i][j]);
                else
                    LHS_p.push_back(LHS[i][j]);
            }

            dynamics::twobody<double> dyn(LHS_p);
            integrator::rk4<double> integrator(&dyn);

            std::vector<double> y_tmp;
            integrator.integrate(tstart,tf,100,LHS_x,y_tmp);
            y.push_back(y_tmp);
        }

        // perform interpolation. For efficiency reason the function that interpolate multiple outputs is used
        // poly will evaluate according to its base
        if(H.size()==0)
            poly.interpolation(LHS,y,H,res_coeffs);
        else{
            poly.solve(H,y,res_coeffs);
        }

        for(int i=0;i<nvar;i++)
            coeffs_all.push_back(res_coeffs[i]);

        x0 = xf;
        tstart=tf;
    }

    end=clock();
    double time = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "Time elapsed : " << time << endl << endl;

    std::ofstream file;
    file.open ("twobody_problem_nonintrusive.txt");

    for(unsigned int k=0; k<coeffs_all.size(); k++){
        for(unsigned int kk=0; kk<coeffs_all[k].size(); kk++)
            file << setprecision(16) << coeffs_all[k][kk] << " ";
        file << "\n";
    }

    file << "\n\n\n\n";
    file.close();

}
