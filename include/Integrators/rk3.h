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
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> x0, std::vector<T> xfinal) const;

        };

    }
}

#endif // RK3_H
