#include "Integrators/base_integrator.h"

using namespace smart;
using namespace integrator;

base_integrator::base_integrator(const std::string &name, const dynamics::base_dynamics *dyn) : m_name(name), m_dyn(dyn){
}


std::string base_integrator::get_name() const{
    return m_name;
}
