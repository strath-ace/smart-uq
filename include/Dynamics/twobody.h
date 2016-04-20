/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTUQ_TWOBODY_H
#define SMARTUQ_TWOBODY_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartuq
{
    namespace dynamics {

        /**
         * @brief The Two-body problem
         *
         * The Two-body problem models the dynamics of an object of mass \f$m\f$ orbiting around Earth is deifned as
         * \f{eqnarray*}{
           \ddot{\mathbf{x}} &=& -\frac{\mu}{r^3}\mathbf{x} + \frac{\mathbf{T}}{m} + \frac{1}{2}\rho \frac{C_D A}{m} \|\mathbf{v}_{rel}\| \mathbf{v}_{rel} + \mathbf{\epsilon}
           \f}
         * where \f$r\f$ is the distance from the Earth,  \f$\mathbf{v}_{rel}\f$ is the Earth relative velocity and the mass of the spacecraft varies as
            \f{eqnarray*}{\dot{m} = - \alpha \|\mathbf{T}\|\f}
         */
        template < class T >
        class twobody: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:

            /**
             * @brief twobody constructor
             *
             * The constructor initialize the problem parameters and scaling factors to the value supplied by the user.
             * The model parameters (10 in total) are supplied in order as:
             * - param 1-3: thrust value alng each direction
             * - param 4: \f$\alpha\f$
             * - param 5-6: \f$(\rho_0,H)\f$ atmospheric parameters (exponential atmospheric model is considered)
             * - param 7: \f$C_DA\f$ produc between aerodynamic coefficient and surface area
             * - param 8-10: unknown acceleration (3 components)
             * Default values for scaling factors are 1 and parameters is a vector of zero values or zero polynomials.
             * @param param vector of parameters (constants value or polynomials)
             * @param t_scale time scaling factor
             * @param r_scale position scaling factor
             */
            twobody(const std::vector<T> &param = std::vector<T>(10), const double &t_scale=1, const double &r_scale=1);

            /**
              * @brief ~twobody deconstructor
              */
            ~twobody();

            /**
             * @brief evaluate evaluate the dinamics of the Two-body problem at a given instant of time and a given state.
             *
             * Function to evaluate the dinamics of the Two-body problem at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            mutable std::vector<T> m_param;
            double m_t_scale;
            double m_r_scale;
            double m_m_scale;

        };

    }
}


#endif // SMARTUQ_TWOBODY_H
