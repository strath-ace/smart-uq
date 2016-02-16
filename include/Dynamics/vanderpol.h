#ifndef VANDERPOL_H
#define VANDERPOL_H

#include "base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace dynamics {

        /**
         * @brief The euler class
         */
        template < class T >
        class vanderpol: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief vanderpol
             * @param param constant force along the y (or z) direction
             */
            vanderpol(const T &mu=.5);

            /**
              * @brief ~vanderpol
              */
            ~vanderpol();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const;

        private:
            T m_mu;
        };

    }
}


#endif // VANDERPOL_H
