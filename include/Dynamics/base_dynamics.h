#ifndef BASE_DYNAMICS_H
#define BASE_DYNAMICS_H

#include <vector>
#include "exception.h"

namespace smart
{
    namespace dynamics {

        /**
         * @brief The base_dynamics class is an abstract class. Any dynamics added to the toolbox needs to inherit from it and implement the method evaluate()
         *
         * The base_dynamics class is an abstract class. Any dynamical system added to the toolbox need to extend this class and implement the method evaluate.
         */
        class base_dynamics
        {

        public:

            /**
             * @brief base_dynamics constructor.
             *
             * In the constructor the name of the dynamics
             * @param name dynamical system name
             */
            base_dynamics(const std::string &name);
            virtual ~base_dynamics();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            virtual int evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const = 0;

            /**
             * @brief get_name
             * @return
             */
            std::string get_name() const;

        protected:
            std::string m_name;

        };

    }
}

#endif // BASE_DYNAMICS_H
