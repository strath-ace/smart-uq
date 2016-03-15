/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef SMART_MATH_H
#define SMART_MATH_H

#include <vector>
#include <cmath>
#include "config.h"


const double ZERO = 1e-15;

template <class T>
T inverse(T x){
    if(fabs(x)<=ZERO){
        std::cout<<"ERROR: Division by zero."<<std::endl;
        throw std::exception();
        //exit(EXIT_FAILURE);
    }
    return 1.0/x;
}

//MATH STUFFS
inline int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline int combination(int n, int k)
{
    int max = std::max(n,k);
    int min = std::min(n,k);
    int res = 1;
    int j = 1;

    while(j<=min){
        res *= max+j;
        j++;
    }

    return res/factorial(min);
}


inline void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, unsigned int count){
    if (count < item.size()){
        for (unsigned int i = 0; i < values.size(); i++) {
            item[count] = values[i];
            unsigned int tmp_count = count + 1;
            rep(res, values, item, tmp_count);
        }
    }else{
        res.push_back(item);
    }
}


inline void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res){
    res.clear();

    std::vector<int> item(k);
    rep(res, values, item, 0);
}

#endif // SMART_MATH_H
