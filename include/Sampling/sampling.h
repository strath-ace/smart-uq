/******************************************************************************
 *                                 SAMPLING_H                                 *
 *                  Sampling Methods of the SMART-UQ toolbox                  *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
----------------- e-mail:  carlos.ortega@strath.ac.uk -----------------------
----------------------- Author:  Carlos Ortega Absil ------------------------
*/
#ifndef SAMPLING_H_
#define SAMPLING_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iomanip>

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
    class sobol
    {
        public:
            sobol(unsigned int dim, unsigned int count = 1);
            ~sobol(){}
            std::vector<T> operator()();
            std::vector<T> operator()(unsigned int n);
        private:
            int i8_bit_lo0 ( long long int n );
            void i8_sobol ( unsigned int dim_num, long long int *seed, T quasi[ ] );
        private:
            unsigned int m_dim; /// Hypercube dimension where sampling
            unsigned int m_count;/// Starting point of the sequence (can be used to skip initial values)
            unsigned int m_dim_num_save;
            bool m_initialized;
            long long int m_maxcol;
            long long int m_seed_save;
            T recipd;
            long long int lastq[1111]; //1111 is maximum dimension.
            long long int poly[1111];
            long long int v[1111][62]; //2^62 is approx. limit of points requested.
    };

    /// Latin Hypercube Sampling
    /**
     * Generates a latin hypersquare sampling
     * in the unit hyper cube [0,1].
     * The code wraps original routines from the link below.
     *
     * NOTE: uses rand and srand. MODIFY WHEN A BETTER RNG IS AVAILABLE
     *
     * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/latin_random/latin_random.html
     * @author carlos.ortega@strath.ac.uk
    */

    template <class T>
    class lhs
    {
        public:
            lhs(unsigned int dim, unsigned int npoints);
            ~lhs(){}
            std::vector<T> operator()();
            std::vector<T> operator()(unsigned int n);
        private:
            std::vector<T> latin_random ( unsigned int dim_num, unsigned int point_num);
            unsigned int *perm_uniform ( unsigned int n);
        private:
            unsigned int m_dim; /// Hypercube dimension where sampling
            unsigned int m_npoints;/// Number of points to generate in set
            bool m_initialised;
            std::vector<T> m_set;
            unsigned int m_next;
    };

}}

#endif /* SAMPLING_H_ */