/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef VANDERPOL_H
#define VANDERPOL_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartuq
{
    namespace dynamics {

        /**
         * @brief The euler class
         */
        template < class T >
        class vanderpol: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief vanderpol
             * @param param constant force along the y (or z) direction
             */
            vanderpol(const T &mu=.5);

            /**
              * @brief ~vanderpol
              */
            ~vanderpol();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            T m_mu;
        };

    }
}


#endif // VANDERPOL_H
