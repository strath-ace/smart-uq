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

#include "../Dynamics/base_dynamics.h"
#include "../exception.h"

namespace smartuq
{
    namespace integrator {

        /**
         * @brief The base_integrator class is a template abstract class. Any fixed stepsize integrator added to the toolbox needs to inherit from it and implement the method integrate()
         *
         * The base_integrator class is a template abstract class. Any fixed stepsize integrator added to the toolbox needs to inherit from it and implement the method that integrates between to given times, initial state and stepsize
         */
        template < class T >
        class base_integrator
        {

        public:

            /**
             * @brief base_integrator constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to eb integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            base_integrator(const std::string &name, const dynamics::base_dynamics<T> *dyn);

            /**
             * @brief ~base_integrator deconstructor
             */
            virtual ~base_integrator();


            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The virtual method is inherited by any subclass. It implements the corresponding integration scheme with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            virtual int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const = 0;

            /**
             * @brief get_name return integrator name
             *
             * Function to get the name of the integration scheme
             * @return
             */
            std::string get_name() const;


        protected:
            /**
             * @brief m_name integrator name
             */
            std::string m_name;
            /**
             * @brief m_dyn pointer to a base_dynamics system to be integrated
             */
            const dynamics::base_dynamics<T> *m_dyn;
        };

    }
}

#endif // BASE_INTEGRATOR_H
