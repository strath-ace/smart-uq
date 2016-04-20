/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_RK3_H
#define SMARTUQ_RK3_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartuq
{
    namespace integrator {

        template <class T> class base_integrator;
        /**
         * @brief The Runge Kutta 3 integrator scheme
         *
         * The class model the Runge Kutta third order integration scheme
         */
        template < class T >
        class rk3: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief rk3 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            rk3(const dynamics::base_dynamics<T> *dyn);

            /**
              * @brief ~rk3 deconstructor
              */
            ~rk3();

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK3 scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const;

        };

    }
}

#endif // SMARTUQ_RK3_H
