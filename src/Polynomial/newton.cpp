/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#include "Polynomial/newton.h"


using namespace smart;
using namespace polynomial;

template < class T >
Newton_Polynomial<T>::Newton_Polynomial(const int &nvar, const int &order) : Polynomial<T>::Polynomial(nvar,order){
    initialize_nodes();
}

template <class T>
Newton_Polynomial<T>::Newton_Polynomial(const int &nvar, const int &order, const int &i) : Polynomial<T>::Polynomial(nvar,order,i){
    initialize_nodes();
}

template <class T>
Newton_Polynomial<T>::Newton_Polynomial(const int &nvar, const int &order, const T &value) : Polynomial<T>::Polynomial(nvar,order,value){
    initialize_nodes();
}

template < class T >
std::string Newton_Polynomial<T>::get_name() const
{
    return "Newton Polynomial";
}

template < class T >
std::string Newton_Polynomial<T>::get_basis_name() const
{
    return "n";
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator+(const Newton_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Newton_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator-(const Newton_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    Newton_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}


//OPERATOR* OVERLOADING FOR DIRECT MULTIPLICATION
template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator*(const Newton_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout<<"NOT IMPLEMENTED"<<std::endl;


    //// MONOMIAL MULTIPLICATION
    // Newton_Polynomial<T> res(m_nvar,m_degree);
    // std::vector<T> res_coeffs(m_coeffs.size());
    // std::vector<T> other_coeffs = other.get_coeffs();
    // int i_0=0;//index offset poly0
    // int j_0=0;//index offset poly1
    // int idx_0=0;//index offset of res
    // for (int deg0=0; deg0<=m_degree; deg0++){ //loop over order of terms of poly0
    //     int max_i=m_J[m_nvar][deg0];
    //     if (deg0>0) i_0=m_N[m_nvar][deg0-1];
    //     for (int deg1=0; deg1<=other.get_degree()-deg0; deg1++){ //loop over orther of terms of poly1
    //         int max_j=m_J[m_nvar][deg1];
    //         if (deg1>0) j_0=m_N[m_nvar][deg1-1];
    //         int deg = deg0+deg1;//order of terms of res
    //         if (deg>0) idx_0=m_N[m_nvar][deg-1];
    //         for (int i=0;i<max_i;i++){
    //             for (int j=0;j<max_j;j++){
    //                 if(fabs(m_coeffs[i_0+i])>ZERO && fabs(other_coeffs[j_0+j])>ZERO){
    //                     //find what term is the result contributing to
    //                     std::vector<int> row0=this->get_row(deg0,i);
    //                     std::vector<int> row1=this->get_row(deg1,j);
    //                     std::vector<int> row(m_nvar);
    //                     for (int k=0;k<m_nvar;k++) row[k]=row0[k]+row1[k];
    //                     int idx = res.get_idx(row);
    //                     //multiply and add contribution
    //                     res_coeffs[idx_0+idx]+= m_coeffs[i_0+i]*other_coeffs[j_0+j];
    //                 }
    //             }
    //         }
    //     }
    // }

    // res.set_coeffs(res_coeffs);
    // return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator/(const Newton_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }



    std::cout<<"NOT IMPLEMENTED"<<std::endl;


}


template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T> Newton_Polynomial<T>::operator-() const{
    
    std::vector<T> coeffs=this->get_coeffs();
    Newton_Polynomial<T> res(m_nvar,m_degree);
    for (int i=0;i<coeffs.size();i++){
        coeffs[i] = -coeffs[i];
    }
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator=(const Newton_Polynomial<T> &other){

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    m_coeffs = other.get_coeffs();
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator+=(const Newton_Polynomial<T> &other){
    *this = Newton_Polynomial<T>::operator+(other);
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator-=(const Newton_Polynomial<T> &other){
    *this = Newton_Polynomial<T>::operator-(other);
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator+=(const T& other){
    *this = Newton_Polynomial<T>::operator+(other);
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator-=(const T& other){
    *this = Newton_Polynomial<T>::operator-(other);
    return *this;
}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator*=(const T& other){
    // *this = Newton_Polynomial<T>::operator*(other);
    // return *this;
    std::cout<<"NOT IMPLEMENTED"<<std::endl;

}

template <class T>
Newton_Polynomial<T>& Newton_Polynomial<T>::operator/=(const T& other){
    // *this = Newton_Polynomial<T>::operator/(other);
    // return *this;
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
}

template <class T>
bool Newton_Polynomial<T>::operator==(const Newton_Polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(m_coeffs==other.get_coeffs())
        return true;
    return false;
}

template <class T>  
bool Newton_Polynomial<T>::operator!=(const Newton_Polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }
    if(m_degree!=other.get_degree()){
        std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(Newton_Polynomial<T>::operator==(other)) return false;
    else return true;

}

//1-d Evaluation method
template <class T>
T Newton_Polynomial<T>::evaluate(const T &x) const {
    if(m_nvar>1){
        std::cout<<"(evaluate) Dimension of point must correspond to number of variables of polynomial."<<std::endl;
        exit(EXIT_FAILURE);
    }
    // if (fabs(x)>1){
    //     std::cout<<"(evaluate) All components of point must belong to [-1,1]."<<std::endl;
    //     exit(EXIT_FAILURE); 
    // }

    return m_coeffs[0]+nested_mult(x,1);

}

//Multivariate Evaluation method
template <class T>
T Newton_Polynomial<T>::evaluate(const std::vector<T> &x) const {
    if(m_nvar!=x.size()){
        std::cout<<"(evaluate) Dimension of point must correspond to number of variables of polynomial."<<std::endl;
        exit(EXIT_FAILURE);
    }
    for (int i=0;i<m_nvar;i++){
        if (fabs(x[i])>1){
            std::cout<<"(evaluate) All components of point must belong to [-1,1]."<<std::endl;
            exit(EXIT_FAILURE); 
        }
    }
    
    std::cout<<"NOT IMPLEMENTED"<<std::endl;
}

//interpolate given a set of values in the nodes (1d)
template <class T>
void Newton_Polynomial<T>::interpolate_nodes(const std::vector<T> &y) {
    if (y.size()!=m_nodes.size()){
        std::cout<<"(interpolate) Values provided must correspond to number of nodes."<<std::endl;
        exit(EXIT_FAILURE);  
    }

    std::vector<T> coeffs = y;
    for (int i = 0 ; i < m_degree ; i++){
        for (int j = m_degree ; j > i ; j--){
            coeffs[j]-=coeffs[j-1];
            coeffs[j]/=(m_nodes[j]-m_nodes[j-i-1]);
        }
    }

    m_coeffs=coeffs;
}


//private routine for evaluation and interpolation
template <class T>
void Newton_Polynomial<T>::initialize_nodes(){
    T pi = 3.141592653589793;
    int n = m_degree;
    m_nodes.resize(0);
    for (int i=0;i<=m_degree;i++){
        m_nodes.push_back(cos(pi*(i+0.5)/(n+1)));
    }
    if (n%2==0) m_nodes[n/2]=0.0;
}

//private routine for 1d evaluation
template <class T>
T Newton_Polynomial<T>::nested_mult(T x, int i) const{
    if (i>m_degree) return 0;
    else return (x-m_nodes[i-1])*(m_coeffs[i]+nested_mult(x,i+1));
}


template class Newton_Polynomial<double>;
template class Newton_Polynomial<float>;
template class Newton_Polynomial<long double>;