/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef BASE_INTEGRATOR_H
#define BASE_INTEGRATOR_H

#include "Dynamics/base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        /**
         * @brief The base_integrator class
         */
        template < class T >
        class base_integrator
        {

        public:

            /**
             * @brief base_integrator
             * @param name
             * @param dyn
             */
            base_integrator(const std::string &name, const dynamics::base_dynamics<T> *dyn);

            /**
             * @brief ~base_integrator
             */
            virtual ~base_integrator();

            /**
             * @brief integrate
             * @param ti
             * @param tend
             * @param nsteps
             * @param x0
             * @param xfinal
             * @return
             */
            virtual int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const = 0;

            /**
             * @brief get_name
             * @return
             */
            std::string get_name() const;


        protected:
            std::string m_name;
            const dynamics::base_dynamics<T> *m_dyn;
        };

    }
}

#endif // BASE_INTEGRATOR_H
