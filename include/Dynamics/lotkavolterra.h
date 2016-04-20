/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_LOTKAVOLTERRA_H
#define SMARTUQ_LOTKAVOLTERRA_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartuq
{
    namespace dynamics {

        /**
         * @brief The Lotka Volterra problem
         *
         * The class implements the 2D differential equation system of the Lotka Volterra problem.
         * Namely
         * \f{eqnarray*}{
            \dot{x} &=& a x - b xy\\
            \dot{y} &=& c y + d xy
          \f}
         * parameters can be set in the constructors.
         */
        template < class T >
        class lotkavolterra: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief lotkavolterra constructor
             *
             * The constructor initializes the parameters of the Lotka-Volterra problem
             * @param[in] param 4 dimensional vector containing the parameters of the problem.
             */
            lotkavolterra(const std::vector<T> &param);

            /**
              * @brief ~lotkavolterra deconstructor
              */
            ~lotkavolterra();

            /**
             * @brief evaluate evaluate the dinamics of the Lotka Volterra problem at a given instant of time and a given state.
             *
             * Function to evaluate the dinamics of the Lotka Volterra problem at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            /**
             * @brief m_param parameters vector
             */
            std::vector<T> m_param;
        };

    }
}

#endif // SMARTUQ_LOTKAVOLTERRA_H
