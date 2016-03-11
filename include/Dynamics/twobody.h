#ifndef AKP_H
#define AKP_H

#include "base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace dynamics {

        /**
         * @brief The twobody class
         */
        template < class T >
        class twobody: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief twobody
             */
            twobody(const std::vector<T> &param = std::vector<T>(10), const double &t_scale=1, const double &r_scale=1);

            /**
              * @brief ~twobody
              */
            ~twobody();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state.
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


#endif // AKP_H
