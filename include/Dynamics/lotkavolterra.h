#ifndef LOTKAVOLTERRA_H
#define LOTKAVOLTERRA_H

#include "base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace dynamics {

        /**
         * @brief The lotka volterra problem class
         */
        template < class T >
        class lotkavolterra: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief lotkavolterra
             * @param param
             */
            lotkavolterra(const std::vector<T> &param=std::vector<T>(4,1));

            /**
              * @brief ~lotkavolterra
              */
            ~lotkavolterra();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            std::vector<T> m_param;
        };

    }
}

#endif // LOTKAVOLTERRA_H
