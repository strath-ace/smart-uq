#include "Dynamics/base_dynamics.h"

using namespace smart;
using namespace dynamics;

base_dynamics::base_dynamics(const std::string &name): m_name(name){
}

std::string base_dynamics::get_name() const{
    return m_name;
}
