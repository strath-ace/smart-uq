/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef RK3_H
#define RK3_H

#include "base_integrator.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 3 integrator class
         */
        template < class T >
        class rk3: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief rk3
             * @param dyn
             */
            rk3(const dynamics::base_dynamics<T> *dyn);

            /**
              * @brief ~rk3
              */
            ~rk3();

            /**
             * @brief integrate
             * @param ti
             * @param tend
             * @param nsteps
             * @param x0
             * @param xfinal
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const;

        };

    }
}

#endif // RK3_H
