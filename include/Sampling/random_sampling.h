/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef RANDOM_SAMPLING_H
#define RANDOM_SAMPLING_H


#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include "base_sampling.h"

using namespace std;

namespace smartuq{
namespace sampling{

    /**
     * @brief Random sampling technique
     *
     * The class implements the random sampling technique
    */

    template <class T>
    class random_sampling : public base_sampling <T>
    {
        private:
            using base_sampling<T>::m_dim;
            using base_sampling<T>::m_name;
            using base_sampling<T>::m_a;
            using base_sampling<T>::m_b;

        public:

        /**
         * @brief random_sampling constructor
         *
         * The constructor initialize the dimension of the hypercube and its ranges
         * @param dim hypercube dimension
         * @param a vector containing the lower bound value of each variable
         * @param b vector containing the upper bound value of each variable
         */
         random_sampling(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b);

         /**
           * @brief ~random deconstructor
           */
         ~random_sampling();

         /**
          * @brief operator () next sample points
          *
          * The method implements the subsequential generation of sample points within the user defined hypercube.
          * Being a random sampling a new random point is generated within the hypercube at each call of the method
          * @return the new sample point
          */
         std::vector<T> operator()() const;

         /**
          * @brief operator () n-th sample point in the sequence
          *
          * Returns the n-th point in the sequence. Being a random sampling
          * a new random point is generated within the hypercube at each call of the method
          * @param[in] n the point along the sequence to be returned
          * @return a vector containing the n-th sample point
          */
         std::vector<T> operator()(const unsigned int &n) const;
    };

}}


#endif // RANDOM_SAMPLING_H
