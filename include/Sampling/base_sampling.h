/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_BASE_SAMPLING_H
#define SMARTUQ_BASE_SAMPLING_H


#include <vector>
#include "../exception.h"

using namespace std;

namespace smartuq{
namespace sampling{


    /**
     * @brief The base_sampling class is a template abstract class. Any sampling added to the toolbox needs to inherit from it and implement the operator ()
     *
     * The base_sampling class is a template abstract class. Any sampling added to the toolbox needs to inherit from it and implement the operator() that generates in sequence the
     * points in the hypercube according to the corresponding scheme.
     */
    template <class T>
    class base_sampling
    {
        public:

            /**
             * @brief base_sampling constructor
             *
             * The constructor initialize the dimension of the hypercube and its span. A name is associate to each sampling technique
             * to be used for error messages.
             * @param dim hypercube dimension
             * @param a vector containing the lower bound value of each variable
             * @param b vector containing the upper bound value of each variable
             * @param name name of the sampling technique
             */
            base_sampling(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b, const string &name);

            /**
             * @brief ~base_sampling deconstructor
             */
            virtual ~base_sampling();

            /**
             * @brief operator () next sample points
             *
             * The method implements the subsequential generation of sample points within the user defined hypercube.
             * It is a virtual method so each class that inherits form it needs to provide an implementation for it.
             * @return the new sample point
             */
            virtual std::vector<T> operator()() const = 0;

            /**
             * @brief operator () n-th sample point in the sequence
             *
             * The method implements the subsequential generation of the n-th sample point within the user defined hypercube.
             * It is a virtual method so each class that inherits form it needs to provide an implementation for it.
             * @param n index of the point to be evaluated
             * @return the n-th sample point according to the scheme
             */
            virtual std::vector<T> operator()(const unsigned int &n) const = 0;

    protected:
            /**
             * @brief map function to map samples from [0,1] to [a,b]
             *
             * @param point sample point to be mapped (must belong to [0,1])
             * @return mapped sample point
             */
            std::vector<T> map(const std::vector<T> &point) const;
    protected:
            /**
             * @brief m_name Sampling technique name
             */
            string m_name;
            /**
             * @brief m_dim Hypercube dimension
             */
            unsigned int m_dim; 
            /**
             * @brief m_a Hypercube range lower bounds
             */
            std::vector<T> m_a;
            /**
             * @brief m_b Hypercube range upper bounds
             */
            std::vector<T> m_b;

    };

}}


#endif // SMARTUQ_BASE_SAMPLING_H
