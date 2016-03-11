/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef BASE_SAMPLING_H
#define BASE_SAMPLING_H


#include <vector>
#include "exception.h"

using namespace std;

namespace smart{
namespace sampling{


    /**
     *@brief base_sampling
     */
    template <class T>
    class base_sampling
    {
        public:

            /**
             * @brief base_sampling
             * @param dim
             * @param name
             */
            base_sampling(const unsigned int &dim, const string &name);

            /**
             * @brief ~base_sampling
             */
            virtual ~base_sampling();

            /**
             * @brief operator ()
             * @return
             */
            virtual std::vector<T> operator()() const = 0;

            /**
             * @brief operator ()
             * @param n
             * @return
             */
            virtual std::vector<T> operator()(const unsigned int &n) const = 0;

        protected:
            string m_name;
            unsigned int m_dim; 

    };

}}


#endif // BASE_SAMPLING_H
