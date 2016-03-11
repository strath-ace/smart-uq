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

namespace smart{
namespace sampling{

    /**
     * @brief
     * @author annalisa.riccardi@strath.ac.uk
    */

    template <class T>
    class random_sampling : public base_sampling <T>
    {
        private:
            using base_sampling<T>::m_dim;
            using base_sampling<T>::m_name;

        public:

            /**
             * @brief random
             * @param dim
             * @param npoints
             */
            random_sampling(const unsigned int &dim);

            /**
              * @brief ~random
              */
            ~random_sampling();

            /**
             * @brief
             *
             * Both return the next point in the sequence (random)
             * @return an std::vector<T> containing the next point
             */
            std::vector<T> operator()() const;
            std::vector<T> operator()(const unsigned int &n) const;
    };

}}


#endif // RANDOM_SAMPLING_H
