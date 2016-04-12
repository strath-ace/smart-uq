/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef BASE_DYNAMICS_H
#define BASE_DYNAMICS_H

#include <vector>
#include "../exception.h"
#include "../Polynomial/smart_polynomial.h"

namespace smartuq
{
    namespace dynamics {

        /**
         * @brief The base_dynamics class is a template abstract class. Any dynamics added to the toolbox needs to inherit from it and implement the method evaluate()
         *
         * The base_dynamics class is an abstract class. Any dynamical system added to the toolbox need to extend this class and implement the method evaluate.
         */
        template < class T >
        class base_dynamics
        {

        public:

            /**
             * @brief base_dynamics constructor.
             *
             * In the constructor the name of the dynamics
             * @param name dynamical system name
             */
            base_dynamics(const std::string &name);

            virtual ~base_dynamics();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            virtual int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const = 0;

            /**
             * @brief get_name
             * @return
             */
            std::string get_name() const;

        protected:
            std::string m_name;

        };

    }
}

#endif // BASE_DYNAMICS_H
