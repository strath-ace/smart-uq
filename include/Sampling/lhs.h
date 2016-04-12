/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef LHS_H
#define LHS_H


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "base_sampling.h"

using namespace std;

namespace smartuq{
namespace sampling{

    /**
     * @brief Latin Hypercube Sampling
     *
     * The class implements the Latin Hypercube sampling technique
     * The code wraps original routines from the link below.
     *
     * @see http://people.sc.fsu.edu/~jburkardt/
    */

    template <class T>
    class lhs : public base_sampling <T>
    {
        private:
            using base_sampling<T>::m_dim;
            using base_sampling<T>::m_name;
            using base_sampling<T>::m_a;
            using base_sampling<T>::m_b;

        public:

            /**
             * @brief lhs constructor
             *
             * The constructor initialize the dimension of the hypercube, its ranges and the number of points requested for the sampling
             * @param dim hypercube dimension
             * @param npoints number of sample points
             * @param a vector containing the lower bound value of each variable
             * @param b vector containing the upper bound value of each variable
             */
            lhs(const unsigned int &dim, const unsigned int &npoints, const std::vector<T>& a, const std::vector<T>& b);

            /**
              * @brief ~lhs deconstructor
              */
            ~lhs();


            /**
             * @brief operator () next sample points
             *
             * The method implements the subsequential generation of sample points within the user defined hypercube.
             * It accesses iteratively the points in the pregenerated latin hypercube
             * @return the new sample point
             */
            std::vector<T> operator()() const;

            /**
             * @brief operator () n-th sample point in the sequence
             *
             * Returns the n-th point in the sequence
             * @param[in] n the point along the sequence to be returned
             * @return a vector containing the n-th sample point
             */
            std::vector<T> operator()(const unsigned int &n) const;

        private:
            std::vector<T> latin_random (const unsigned int &dim_num, const unsigned int &point_num) const;
            unsigned int *perm_uniform (const unsigned int &n) const;

        private:
            unsigned int m_npoints;
            mutable bool m_initialised;
            mutable std::vector<T> m_set;
            mutable unsigned int m_next;
    };

}}


#endif // LHS_H
