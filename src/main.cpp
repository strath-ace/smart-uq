/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include <iostream>

#include "../include/smartuq.h"

using namespace std;

int main(){

    cout << "Welcome to SMART-UQ!" << endl;

    int nvar = 3;
    int ndegree = 15;

    taylor_polynomial<double> x1(nvar,ndegree,0, 1.0, 2.0);
    taylor_polynomial<double> x2(nvar,ndegree,1, 1.0, 2.0);
    taylor_polynomial<double> x3(nvar,ndegree,2, 1.0, 2.0);

    taylor_polynomial<double> f = x1*x1 + x2 + x1*x2*x3;

    taylor_polynomial<double> res1 = sqrt(f*f) - f;
    taylor_polynomial<double> res2 = exp(log(f)) - f;
    taylor_polynomial<double> res3 = log(exp(f)) - f;
    taylor_polynomial<double> res4 = sin(f)*sin(f) + cos(f)*cos(f) - 1.0;
    taylor_polynomial<double> res5 = f/f -1.0;

    cout << "Welcome to SMART-UQ!" << endl;
    return 0;

}
