#ifndef EULER_H
#define EULER_H

#include "base_integrator.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        /**
         * @brief The euler class
         */
        class euler: public base_integrator
        {

        public:

            /**
             * @brief euler
             * @param dyn
             */
            euler(const dynamics::base_dynamics *dyn);

            /**
              * @brief ~euler
              */
            ~euler();

            /**
             * @brief integrate
             * @param ti
             * @param tend
             * @param nsteps
             * @param x0
             * @param xfinal
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<double> x0, std::vector<double> xfinal) const;
	
        };

    }
}

#endif // EULER_H
