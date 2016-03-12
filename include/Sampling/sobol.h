/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SOBOL_H
#define SOBOL_H

#include "constants.h"
#include "base_sampling.h"
#include <sstream>

using namespace std;

namespace smart{
namespace sampling{

    /**
     * @brief Sobol quasi-random point sequence
     *
     * Generates a quasi-random sequence of points
     * in the predefined hypercube using the Sobol sequence.
     * The code wraps original routines from the link below.
     *
     * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html
    */
    template <class T>
    class sobol : public base_sampling <T>
    {
        private:
            using base_sampling<T>::m_dim;
            using base_sampling<T>::m_name;
            using base_sampling<T>::m_a;
            using base_sampling<T>::m_b;

        public:

            /**
            * @brief sobol constructor
            *
            * The constructor initialize the dimension of the hypercube, its ranges and the starting point for sequence generation
            * @param dim dimension of the hypercube
            * @param a vector containing the lower bound value of each variable
            * @param b vector containing the upper bound value of each variable
            * @param count starting point of the sequence. choosing 0 will add the point x=0 (initial point need to be given in the unitary hypercube [0,1])
            * @throws value_error if dim not in [1,1111]
            */
            sobol(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b, const unsigned int &count = 1);

            /**
              * @brief ~sobol deconstructor
              */
            ~sobol();

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
            int i8_bit_lo0 ( long long int n ) const;
            void i8_sobol ( unsigned int dim_num,  long long int *seed,  T quasi[ ] ) const;

        private:
            mutable unsigned int m_count;/// Starting point of the sequence (can be used to skip initial values)
            mutable unsigned int m_dim_num_save;
            mutable bool m_initialized;
            mutable long long int m_maxcol;
            mutable long long int m_seed_save;
            mutable T recipd;
            mutable long long int lastq[1111]; //1111 is maximum dimension.
            mutable long long int poly[1111];
            mutable long long int v[1111][62]; //2^62 is approx. limit of points requested.
    };

}}


#endif // SOBOL_H
