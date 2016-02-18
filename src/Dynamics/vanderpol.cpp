#include "Dynamics/vanderpol.h"


using namespace smart;
using namespace dynamics;

template < class T >
vanderpol<T>::vanderpol(const T &mu) : base_dynamics<T>("Van der Pol dynamical system"), m_mu(mu)
{

}

template < class T >
vanderpol<T>::~vanderpol()
{

}

template < class T >
int vanderpol<T>::evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_exception(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=2)
        smart_exception(m_name+": the state dimension needs to be 2");
    if(state.size()!=dstate.size())
        smart_exception(m_name+": state and dstate must have the same dimension");

    dstate[0] = state[1];
    dstate[1] = m_mu*(1.0-state[0]*state[0])*state[1]-state[0];

    return 0;
}




