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

using namespace smart;

int main(){

    /**
      17 uncertain variables (7 states + 10 parameters)
      4 degree polynomial expansion
      **/

    //EXAMPLE INPUT
    bool scale_problem = true; //true for non-dimensional problem
    bool print_results_to_file = true;
    bool print_time_to_screen = true;

    //algebra
    int nvar = 7;
    int nparam  = 10;
    int poly_degree = 4;

    //polynomial allocation
    std::vector<taylor_polynomial<double> > x0, param, xf;
    
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

    //initialise 7 state variables as Taylor base of order 1 in the variable i mapped to [x-unc_x, x+unc_x]
    for(int i=0;i<nvar;i++){
        x0.push_back(taylor_polynomial<double>(nvar+nparam, poly_degree, i, x[i]-unc_x[i], x[i]+unc_x[i]));
    }
    //initialise 10 parameter variables as Taylor base of order 1 in the variable 7+i mapped to [p-unc_p, p+unc_p]
    for(int i=0;i<nparam;i++){
        param.push_back(taylor_polynomial<double>(nvar+nparam, poly_degree, 7+i, p[i]-unc_p[i], p[i]+unc_p[i]));
    }

    //timing
    clock_t begin, end;
    begin=clock();
    
    //fast multiplication
    x0[0].initialize_M(nvar+nparam,poly_degree);    

    //dynamical system
    dynamics::twobody < taylor_polynomial<double> > dyn(param, t_scale, r_scale);
    integrator::rk4 < taylor_polynomial<double> > integrator(&dyn);

    //propagation (MAIN LOOP)
    std::vector<std::vector<double> > coeffs_all;
    deltat = 1000.0 / t_scale; //snapshot frequency
    for(int i=0; i+1<tend/deltat; i++){
        tf += deltat;
        integrator.integrate(tstart,tf,100,x0,xf);
        x0 = xf;
        tstart=tf;
        //save expansions
        for(int j=0; j<nvar; j++){
            std::vector<double> coeffs = xf[j].get_coeffs();
            coeffs_all.push_back(coeffs);
        }
    }

    //deallocation of heavy structures
    x0[0].delete_M();

    //timing
    end=clock();
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    if (print_time_to_screen) cout << "example_taylor, time elapsed : " << time_elapsed << endl;

    //printing
    if(print_results_to_file){
        std::ofstream file;
        file.open ("twobody_problem_taylor.txt");
        for(unsigned int k=0; k<coeffs_all.size(); k++){
            for(unsigned int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
        }
        file << "\n\n\n\n";
        file.close();
    }
}
