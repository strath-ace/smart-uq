#include "Dynamics/twobody.h"


using namespace smart;
using namespace dynamics;

template < class T >
twobody<T>::twobody(const std::vector<T> &param, const double &t_scale, const double &r_scale, const double &m_scale) : base_dynamics<T>("Two Body Problem"),
    m_param(param), m_t_scale(t_scale), m_r_scale(r_scale), m_m_scale(m_scale)
{
    if(m_param.size()!=10)
        smart_exception(m_name+": the parameters list need to be of size 10");

}

template < class T >
twobody<T>::~twobody()
{

}

template < class T >
int twobody<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_exception(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=7)
        smart_exception(m_name+": the state dimension needs to be 7");

    std::vector<T> sstate(state);
    dstate.clear();

    //scaling of states
    for(int i=0;i<3;i++){
        sstate[i] /= m_r_scale;
        sstate[i+3] /= m_r_scale/m_t_scale;
    }
    sstate[6] /= m_m_scale;
    //scaling of parameters
    m_param[0] /= m_m_scale*m_r_scale/(m_t_scale*m_t_scale); //thrust
    m_param[1] /= m_m_scale*m_r_scale/(m_t_scale*m_t_scale); //thrust
    m_param[2] /= m_m_scale*m_r_scale/(m_t_scale*m_t_scale); //thrust
    m_param[3] /= m_t_scale/m_r_scale; //alpha
    m_param[4] /= m_m_scale/pow(m_r_scale,3); //rho
    m_param[5] /= m_r_scale; //H
    m_param[6] /= pow(m_r_scale,2); //CDA
    m_param[7] /= m_r_scale/pow(m_t_scale,2); //epislon
    m_param[8] /= m_r_scale/pow(m_t_scale,2); //epislon
    m_param[9] /= m_r_scale/pow(m_t_scale,2); //epislon

    double radius_earth = 6378*pow(10,3)/m_r_scale;
    double mu = 398600.4415*pow(10,9)/( pow(m_r_scale,3)/pow(m_t_scale,2) );

    T r = sqrt(sstate[0]*sstate[0]+sstate[1]*sstate[1]+sstate[2]*sstate[2]);
    T tmp_3D =  mu/pow(r, 3);

    // atmospheric model
    T rho = m_param[4]*exp(-(r-radius_earth-900000/m_r_scale)/m_param[5]);

    //relative velocity
    T rel_v_x = sstate[3]-7.2921150*m_t_scale*pow(10,-5)*sstate[1];
    T rel_v_y = sstate[4]+7.2921150*m_t_scale*pow(10,-5)*sstate[0];
    T mod_rel_v = sqrt(rel_v_x*rel_v_x+rel_v_y*rel_v_y+sstate[5]*sstate[5]);

    //drag computation
    T tmp_drag = 0.5*rho*m_param[6]*mod_rel_v/sstate[6];

    dstate.push_back(sstate[3]);
    dstate.push_back(sstate[4]);
    dstate.push_back(sstate[5]);
    dstate.push_back(-tmp_3D*sstate[0]+m_param[0]/sstate[6]+m_param[7]-tmp_drag*rel_v_x);
    dstate.push_back(-tmp_3D*sstate[1]+m_param[1]/sstate[6]+m_param[8]-tmp_drag*rel_v_y);
    dstate.push_back(-tmp_3D*sstate[2]+m_param[2]/sstate[6]+m_param[9]-tmp_drag*sstate[5]);
    dstate.push_back(-m_param[3]*sqrt(m_param[0]*m_param[0]+m_param[1]*m_param[1]+m_param[2]*m_param[2]));

    return 0;

}


template class twobody<double>;
template class twobody<float>;
template class twobody<long double>;
template class twobody<polynomial::chebyshev_polynomial<double> >;
template class twobody<polynomial::chebyshev_polynomial<float> >;
template class twobody<polynomial::chebyshev_polynomial<long double> >;
template class twobody<polynomial::taylor_polynomial<double> >;
template class twobody<polynomial::taylor_polynomial<float> >;
template class twobody<polynomial::taylor_polynomial<long double> >;

