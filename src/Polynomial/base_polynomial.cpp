/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#include "Polynomial/base_polynomial.h"

using namespace smart;
using namespace polynomial;



/**************/
/*CONSTRUCTORS*/
/**************/
template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_manipulated_to_monomial(false), m_J(0), m_N(0){

    if(vars<0){
        smart_exception(m_name+"Polynomials need to have a positive number of variables");
    }
    if(order<0){
        smart_exception(m_name+"Polynomials need to have a positive order");
    }

    int n = combination(vars,order);

    if(n>smart::constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_exception(m_name+"The size of the algebra is too big. Reduce polynomial order rnumber of variables. Youcan incur in memory issues");
    }

    m_coeffs.resize(n);

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();

}

template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order, const int &i): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_manipulated_to_monomial(false), m_J(0), m_N(0){

    if(vars<0){
        smart_exception(m_name+"Polynomials need to have a positive number of variables");
    }
    if(order<=0){
        smart_exception(m_name+"Polynomials need to have a positive order");
    }
    if(i<0 || i>=vars){
        smart_exception(m_name+"First order Polynomial constructor need a variable index between [0,nvars-1]");
    }

    //allocate memory for coefficients vector

    int n = combination(vars,order);

    if(n>smart::constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_exception(m_name+"The size of the algebra is too big. Reduce polynomial order rnumber of variables. Youcan incur in memory issues");
    }

    m_coeffs.resize(n);
    m_coeffs[i+1] = 1.0;

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();

}


template < class T >
base_polynomial<T>::base_polynomial(const int &vars, const int &order, const T &value): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_manipulated_to_monomial(false), m_J(0), m_N(0){

    if(vars<0){
        smart_exception(m_name+"Polynomials need to have a positive number of variables");
    }
    if(order<=0){
        smart_exception(m_name+"Polynomials need to have a positive order");
    }

    //allocate memory for coefficients vector

    int n = combination(vars,order);

    if(n>smart::constants::MAX_POLYNOMIAL_ALGEBRA_SIZE){
        smart_exception(m_name+"The size of the algebra is too big. Reduce polynomial order rnumber of variables. Youcan incur in memory issues");
    }

    m_coeffs.resize(n);
    m_coeffs[0] = value;

    //save some info
    m_degree = order;
    m_nvar = vars;

    m_J.resize(vars+1);
    m_N.resize(vars+1);
    for(int i=0; i<=vars; i++){
        m_J[i].resize(order+1);
        m_N[i].resize(order+1);
    }

    initialize_J();
    initialize_N();


}

template <class T>
void base_polynomial<T>::interpolation(const std::vector<std::vector<T> > &x, const std::vector<T> &y) const{
    if(x.size()==0)
        smart_exception(m_name+"for polynomial interpolation non empty nodal values need to be provided");
    if(x.size()!=y.size())
        smart_exception(m_name+"for polynomial interpolation, the number of nodes and nodal values need to be the same");
    if(x[0].size()!=m_nvar)
        smart_exception(m_name+"the number of variables is not the same as in the algebra");

    int npoints = x.size();
    int ncoeffs = m_coeffs.size();

    if(npoints<ncoeffs)
        smart_exception(m_name+"the number of interpolation point need to be equal or greater than the size of the algebra");

    Eigen::MatrixXd base_matrix (npoints,ncoeffs);
}



/******************************/
/*Monomial multiplication     */
/******************************/
//initialisation of static terms for multiplication
template <class T>
std::vector<int> base_polynomial<T>::m_M(1,0);
template <class T>
int base_polynomial<T>::m_Mnvar = 0;
template <class T>
int base_polynomial<T>::m_Mdegree = 0;
template <class T>
void base_polynomial<T>::monomial_multiplication(const base_polynomial<T> &x1, const base_polynomial<T> &x2, base_polynomial<T> &res_poly) const{

    if(!x1.m_manipulated_to_monomial || !x2.m_manipulated_to_monomial){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(x1.get_nvar()!=x2.get_nvar() || x1.get_nvar()!=res_poly.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(x1.get_degree()!=x2.get_degree() || x1.get_degree()!=res_poly.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int degree = x1.get_degree();
    int nvar = x1.get_nvar();
    std::vector<std::vector<int> > J = x1.get_J();
    std::vector<std::vector<int> > N = x1.get_N();

    // use M instead of searching index for faster multiplication
    bool use_M = false;
    if(nvar == m_Mnvar && degree == m_Mdegree) use_M = true;

    std::vector<T> x1_coeffs = x1.get_coeffs();
    std::vector<T> x2_coeffs = x2.get_coeffs();
    std::vector<T> res_coeffs(x1_coeffs.size());
    int i_0, j_0, idx_0; //index offsets
    int idx;
    int coeff=0;


    for (int deg0=0; deg0<=degree; deg0++){ //loop over order of terms of poly0
        int max_i=J[nvar][deg0];
        if (deg0==0) i_0=0;
        else i_0=N[nvar][deg0-1];

        for (int deg1=0; deg1<=degree-deg0; deg1++){ //loop over other of terms of poly1
            int max_j=J[nvar][deg1];
            if (deg1==0) j_0=0;
            else j_0=N[nvar][deg1-1];

            int deg = deg0+deg1;//order of terms of res
            if (deg==0) idx_0=0;
            else idx_0=N[nvar][deg-1];

            for (int i=0;i<max_i;i++){
                for (int j=0;j<max_j;j++){
                    if(fabs(x1_coeffs[i_0+i])>ZERO && fabs(x2_coeffs[j_0+j])>ZERO){
                        if (use_M) idx = m_M[coeff];
                        else{
                        //find what term is the result contributing to
                            std::vector<int> row0=x1.get_row(i,deg0);
                            std::vector<int> row1=x1.get_row(j,deg1);
                            std::vector<int> row(nvar);
                            for (int k=0;k<nvar;k++) row[k]=row0[k]+row1[k];
                            idx = res_poly.get_idx(row);
                        }
                        //multiply and add contribution
                        res_coeffs[idx_0+idx]+= x1_coeffs[i_0+i]*x2_coeffs[j_0+j];
                    }
                    coeff++;
                }
            }

        }
    }

    res_poly.set_coeffs(res_coeffs);

}


/******************************/
/*J & N matrix manipulation   */
/******************************/
template < class T >
void base_polynomial<T>::initialize_J(){
    int i,j;

    //fill J
    for(j = 0; j <= m_degree; ++j)
        m_J[1][j] = 1;
    for(i = 2; i <= m_nvar; ++i) {
        m_J[i][0] = 1;
        for(j = 1; j <= m_degree; ++j)
            m_J[i][j] = m_J[i][j-1]+m_J[i-1][j];
    }

}

template < class T >
void base_polynomial<T>::initialize_N(){
    int i,j;
    //fill N
    for(i = 1; i <= m_nvar; ++i)
        m_N[i][0] = m_J[i][0];
    for(i = 1; i <= m_nvar; ++i) {
        for(j = 1; j <= m_degree; ++j)
            m_N[i][j] = m_N[i][j-1]+m_J[i][j];
    }

}

template < class T >
std::vector<int> base_polynomial<T>::get_row(const int &idx, const int &deg) const{

    int i, j, l;
    int idx_tmp = idx;
    std::vector<int> k(m_nvar);

    if(deg > 0)
        idx_tmp += m_N[m_nvar][deg-1];
    l = m_nvar;
    for(i = 0; i < m_nvar; ++i) {
        if(idx_tmp == 0) {
            for(j = i; j < m_nvar; k[j++] = 0);
            return k;
        }
        if(i == 0)
            j = deg;
        else {
            for(j = 0; idx_tmp >= m_N[l][j]; ++j);
            k[i-1] -= j;
        }
        k[i] = j;
        idx_tmp -= m_N[l--][j-1];
    }

    return k;
}

template < class T >
int base_polynomial<T>::get_idx(const std::vector<int> &k) const{
    int i, j, l;

    int idx = 0;
    for(i = m_nvar-1; (i>=0) && (k[i]==0); --i);

    if(i < 0)
        return idx;
    j = m_nvar-i;
    l = -1;
    while(i >= 0) {
        l += k[i--];
        idx += m_N[j++][l];
    }
    idx -= m_N[m_nvar][l];

    return idx;

}

//initialisation of M
template <class T>
void base_polynomial<T>::initialize_M(const int &nvar, const int &degree){

    base_polynomial<T> poly(nvar,degree);
    std::vector<int> J=poly.get_J()[nvar];
    std::vector<int> N=poly.get_N()[nvar];
    std::vector<int> M;
    for (int deg0=0; deg0<=degree; deg0++){ //loop over order of terms of poly0
        int max_i=J[deg0];
        for (int deg1=0; deg1<=degree-deg0; deg1++){ //loop over other of terms of poly1
            int max_j=J[deg1];
            //int deg = deg0+deg1; //order of terms of result
            for (int i=0;i<max_i;i++){
                for (int j=0;j<max_j;j++){
                    //find what term is the result contributing to
                    std::vector<int> row0=poly.get_row(i,deg0);
                    std::vector<int> row1=poly.get_row(j,deg1);
                    std::vector<int> row(nvar);
                    for (int k=0;k<nvar;k++) row[k]=row0[k]+row1[k];
                    //store it to avoid computing it in every multiplication
                    M.push_back(poly.get_idx(row));
                }
            }
        }
    }
    //modify static stuff
    m_M=M;
    m_Mdegree = degree;
    m_Mnvar = nvar;
}

template <class T>
void base_polynomial<T>::delete_M(){
    m_M.resize(1);
    m_Mnvar=0;
    m_Mdegree=0;
}


template class base_polynomial<double>;
template class base_polynomial<float>;
template class base_polynomial<long double>;
