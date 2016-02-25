#include "Dynamics/lotkavolterra.h"


using namespace smart;
using namespace dynamics;

template < class T >
lotkavolterra<T>::lotkavolterra(const std::vector<T> &param) : base_dynamics<T>("Lotka-Volterra dynamical system"), m_param(param)
{
    //sanity checks
    if(m_param.size()!=4)
        smart_exception(m_name+": the size of the parameters vector need to be 4");

}

template < class T >
lotkavolterra<T>::~lotkavolterra()
{

}

template < class T >
int lotkavolterra<T>::evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_exception(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=2)
        smart_exception(m_name+": the state dimension needs to be 2");
    if(state.size()!=dstate.size())
        smart_exception(m_name+": state and dstate must have the same dimension");

    T xy = state[0]*state[1];
    dstate[0] = m_param[0]*state[0]-m_param[1]*xy;
    dstate[1] = -m_param[2]*state[1]+m_param[3]*xy;

    return 0;
}



