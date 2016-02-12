#ifndef BASE_INTEGRATOR_H
#define BASE_INTEGRATOR_H

#include "Dynamics/base_dynamics.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        class base_integrator
        {

        public:


            base_integrator(const std::string &name, const dynamics::base_dynamics *dyn);
            virtual ~base_integrator();

            virtual int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<double> x0, std::vector<double> xfinal) const = 0;

            std::string get_name() const;


        protected:
            std::string m_name;
	    const dynamics::base_dynamics *m_dyn;	
        };

    }
}

#endif // BASE_INTEGRATOR_H
