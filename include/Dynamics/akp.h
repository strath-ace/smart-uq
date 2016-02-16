#ifndef AKP_H
#define AKP_H

#include "base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace dynamics {

        /**
         * @brief The euler class
         */
        template < class T >
        class akp: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief akp
             * @param dim dimension of the problem (2 and 3 allowed)
             * @param force constant force along each direction
             */
            akp(const int &dim=2, const std::vector<T> &force=0);

            /**
              * @brief ~akp
              */
            ~akp();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const;

        private:
            int m_dim;
            std::vector<T> m_force;
        };

    }
}


#endif // AKP_H
