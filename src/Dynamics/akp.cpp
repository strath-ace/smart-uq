#include "Dynamics/akp.h"


using namespace smart;
using namespace dynamics;

template < class T >
akp<T>::akp(const int &dim=2, const std::vector<T> &force=0) : base_dynamics("Accelerated Kepler Problem"), m_dim(dim), m_force(force)
{
    if(force.size()==0)
        m_force = std::vector<T>(dim);

    //sanity checks
    if(m_dim!=2 || m_dim!=3)
        smart_exception(m_name+": 2D and 3D Kepler problem supported");
    if(m_force.size()!=m_dim)
        smart_exception(m_name+": the force vector when supplied must have the same size as the dimension of the problem (one force along each direction)");

}

template < class T >
akp<T>::~akp()
{

}

template < class T >
int akp<T>::evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_exception(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=m_dim)
        smart_exception(m_name+": the state dimension needs to be the same as problem dimension");
    if(state.size()!=dstate.size())
        smart_exception(m_name+": state and dstate must have the same dimension");


    T mu = 1.0;
    T norm2 = 0;

    //computing square norm of distance (computed once and for all to save computations)
    for(int i=0; i<m_dim ; i++){
        norm2 += state[i]*state[i];
    }
    T tmp =  mu/pow(sqrt(norm2), 3);

    //evaluating system
    for(int i=0; i<m_dim ; i++){
        dstate[i] = state[i];
        dstate[i+m_dim] = -1.0*tmp*state[i]+m_force[i];
    }

    return 0;

}


