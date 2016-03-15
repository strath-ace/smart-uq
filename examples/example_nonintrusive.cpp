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

    //EXAMPLE INPUT
    bool monomial_base = true; //true for output coefficients in monomial base, false for chebyshev
    bool scale_problem = true; //true for non-dimensional problem
    bool print_results_to_file = true;
    bool print_time_to_screen = true;

    //algebra
    int nvar = 7;
    int nparam  = 10;
    int poly_degree = 4;

    //number of points in the sample
    int nsamples = combination(nvar+nparam,poly_degree);

    //scaling of fundamental units
    double m_scale = scale_problem ? 2000.0 : 1.0; //M0_spacecraft
    double r_scale = scale_problem ? 6378136.0 : 1.0; //DU_Earth
    double t_scale = scale_problem ? 806.78 : 1.0; //TU_Earth

    //initialisation problem constants
    double sma = 7378.0*pow(10,3) / r_scale;
    double period = 2.0*M_PI/pow(sma * r_scale,-3.0/2.0)/sqrt(398600.4415*pow(10,9)) / t_scale;
    double tstart = 0.0;
    double tf = 0.0;
    double deltat = 0.0;
    double tend = period;

    //initialisation: allocations
    std::vector<double> x(nvar), p(nparam), unc_x(nvar), unc_p(nparam);

    //initialisation: nominal initial states
    x[0] = 7338.0*pow(10,3) / r_scale; //x
    x[1] = 0.0; //y
    x[2] = 0.0; //z
    x[3] = 0.0; //v_x
    x[4] = 2.0*M_PI*sma/period; //v_y
    x[5] = 0.0; //v_z
    x[6] = 2000.0 / m_scale; //m

    //initialisation: nominal parameters
    p[0] = 0.0; //T_x
    p[1] = 0.5 / (m_scale*r_scale / pow(t_scale,2)); //T_y
    p[2] = 0.0; //T_z
    p[3] = 1.0/30000.0 / (t_scale / r_scale); //specific fuel consumption
    p[4] = 5.245*pow(10,-15) / (m_scale/pow(r_scale,3)); //rho_0
    p[5] = 181050.0 / r_scale; //H
    p[6] = 4.4 / pow(r_scale,2); //C_D*A
    p[7] = 0.0; //epsilon_x
    p[8] = 0.0; //epsilon_y
    p[9] = 0.0; //epsilon_z

    //initialisation: uncertainty in initial states
    unc_x[0] = 1000.0 / r_scale; 
    unc_x[1] = 1000.0 / r_scale;
    unc_x[2] = 1000.0 / r_scale;
    unc_x[3] = 5.0 / (r_scale/t_scale);
    unc_x[4] = 5.0 / (r_scale/t_scale);
    unc_x[5] = 5.0 / (r_scale/t_scale);
    unc_x[6] = 1.0 / m_scale;

    //initialisation: uncertainty in parameters
    unc_p[0] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[1] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[2] = 0.05*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    unc_p[3] = 0.05*p[3];
    unc_p[4] = 0.01*p[4];
    unc_p[5] = 0.01*p[5];
    unc_p[6] = 0.01*p[6];
    unc_p[7] = 0.0001 / (r_scale/pow(t_scale,2));
    unc_p[8] = 0.0001 / (r_scale/pow(t_scale,2));
    unc_p[9] = 0.0001 / (r_scale/pow(t_scale,2));

    //ranges for sampling
    std::vector<double> ranges_lb(nvar+nparam), ranges_ub(nvar+nparam);
    for(int i=0; i<nvar; i++){
        ranges_lb[i] = x[i]-unc_x[i];
        ranges_ub[i] = x[i]+unc_x[i];
    }
    for(int i=0; i<nparam; i++){
        ranges_lb[nvar+i] = p[i]-unc_p[i];
        ranges_ub[nvar+i] = p[i]+unc_p[i];
    }

    //timing
    clock_t begin,end;
    begin=clock();

    //construct LHS sampling
    sampling::lhs<double> lhs_gen(nvar+nparam,nsamples,ranges_lb, ranges_ub);
    std::vector<std::vector<double> > LHS, H;
    for(int i=0;i<nsamples;i++){
        std::vector<double> sample=lhs_gen();
        LHS.push_back(sample);
    }

    //initialise polynomial for interpolation. the monomial flag will define the base
    chebyshev_polynomial<double> poly(nvar+nparam,poly_degree, ranges_lb, ranges_ub, monomial_base);

    //propagation (MAIN LOOP)
    std::vector<std::vector<double> > coeffs_all;
    deltat = 1000.0 / t_scale;
    std::vector<std::vector<double> > y;

    for(int t=0; t+1<tend/deltat; t++){
        std::vector<std::vector<double> > res_coeffs;
        tf += deltat;

        // perform nsamples forward integrations
        for(int i=0;i<nsamples;i++){
            std::vector<double> LHS_p, LHS_x;
            // separate states from parameters in the LHS sampling
            for(int j=0;j<nvar+nparam; j++){
                if(j<nvar){
                    if(t==0)
                       LHS_x.push_back(LHS[i][j]);
                }
                else{
                    LHS_p.push_back(LHS[i][j]);
                }
            }
            if(t>0)
                LHS_x = y[i];

            //initialise dynamics and integrator with corresponding set of parameters
            dynamics::twobody<double> dyn(LHS_p, t_scale, r_scale);
            integrator::rk4<double> integrator(&dyn);

            //perform integration and save results
            std::vector<double> y_tmp;
            integrator.integrate(tstart,tf,100,LHS_x,y_tmp);
            if(t==0)
                y.push_back(y_tmp);
            else
                y[i] = y_tmp;

        }

        // perform interpolation - for efficiency reasons the function that interpolates multiple outputs is used
        if(H.size()==0)
            poly.interpolation(LHS,y,H,res_coeffs);
        else{
            poly.solve(H,y,res_coeffs);
        }

        for(int i=0;i<nvar;i++)
            coeffs_all.push_back(res_coeffs[i]);

        tstart=tf;
    }

    //timing
    end=clock();
    double time = (double (end-begin))/CLOCKS_PER_SEC;
    if(print_time_to_screen) cout << "example_nonintrusive, time elapsed : " << time << endl << endl;

    //printing
    if(print_results_to_file){
        std::ofstream file;
        file.open ("twobody_problem_nonintrusive_monomial.txt");
        for(unsigned int k=0; k<coeffs_all.size(); k++){
            for(unsigned int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
        }
        file << "\n\n\n\n";
        file.close();
    }

}
