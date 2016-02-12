#ifndef SOBOL_H
#define SOBOL_H

/******************************************************************************
 *                                 SOBOL_H                                    *
 *                  Sobol Sampling Method of the SMART-UQ toolbox             *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
----------------- e-mail: carlos.ortega@strath.ac.uk ------------------------
----------------------- Author: Carlos Ortega Absil -------------------------
*/

#include "initializers.h"
#include "base_sampling.h"
#include <sstream>

using namespace std;

namespace smart{
namespace sampling{

    /// Sobol quasi-random point sequence
    /**
     * Generates a quasi-random sequence of points
     * in the unit hyper cube [0,1] using the Sobol sequence.
     * The code wraps original routines from the link below.
     *
     * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html
     * @author carlos.ortega@strath.ac.uk
    */
    template <class T>
    class sobol : public base_sampling <T>
    {
        private:
            using base_sampling<T>::m_dim;
            using base_sampling<T>::m_name;

        public:

            /**
            * @param[in] dim dimension of the hypercube
            * @param[in] count starting point of the sequence. choosing 0 wil add the point x=0
            * @throws value_error if dim not in [1,1111]
            */

            sobol(const unsigned int &dim, const unsigned int &count = 1);


            ~sobol();

            /**
             * Returns the next point in the sequence
             *
             * @return an std::vector<T> containing the next point
             */

            std::vector<T> operator()() const;

            /**
             * Returns the n-th point in the sequence
             *
             * @param[in] n the point along the sequence to be returned
             * @return an std::vector<T> containing the n-th point
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
