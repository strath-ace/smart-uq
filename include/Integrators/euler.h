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
        template < class T >
        class euler: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief euler
             * @param dyn
             */
            euler(const dynamics::base_dynamics<T> *dyn);

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
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> x0, std::vector<T> xfinal) const;
	
        };

    }
}

#endif // EULER_H
