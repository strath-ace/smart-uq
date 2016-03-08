#ifndef RK4_H
#define RK4_H

#include "base_integrator.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 4 integrator class
         */
        template < class T >
        class rk4: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief rk4
             * @param dyn
             */
            rk4(const dynamics::base_dynamics<T> *dyn);

            /**
              * @brief ~rk4
              */
            ~rk4();

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

#endif // RK4_H
