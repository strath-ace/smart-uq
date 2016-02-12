#ifndef EULER_H
#define EULER_H

#include "base_integrator.h"
#include "exception.h"

namespace smart
{
    namespace integrator {

        class euler: public base_integrator
        {

        public:

            euler(const dynamics::base_dynamics *dyn);
            ~euler();

            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<double> x0, std::vector<double> xfinal) const;
	
        };

    }
}

#endif // EULER_H
