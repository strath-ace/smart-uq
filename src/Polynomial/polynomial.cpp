/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#include "Polynomial/polynomial.h"

using namespace smart;
using namespace polynomial_algebra;



/**************/
/*CONSTRUCTORS*/
/**************/
template < class T >
polynomial<T>::polynomial(const int &vars, const int &order): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_J(0), m_N(0){

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
polynomial<T>::polynomial(const int &vars, const int &order, const int &i): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_J(0), m_N(0){

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
polynomial<T>::polynomial(const int &vars, const int &order, const T &value): m_name("Polynomial"), m_coeffs(0), m_degree(0), m_nvar(0),
    m_J(0), m_N(0){

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


/******************************/
/*ARITHMETIC OPERATIONS (+-*) */
/******************************/
template < class T >
polynomial<T> polynomial<T>::operator+(const polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }

    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;


}
template < class T >
polynomial<T> polynomial<T>::operator-(const polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

//OPERATOR* OVERLOADING FOR DIRECT MULTIPLICATION
//initialisation of static terms for multiplication
template <class T>
std::vector<int> polynomial<T>::m_M(1,0);
template <class T>
int polynomial<T>::m_Mnvar = 0;
template <class T>
int polynomial<T>::m_Mdegree = 0;
template <class T>
polynomial<T> polynomial<T>::operator*(const polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    // use M instead of searching index for faster multiplication
    bool use_M = false;
    if(m_nvar == m_Mnvar && m_degree == m_Mdegree) use_M = true;

    polynomial<T> res(m_nvar,m_degree);
    std::vector<T> res_coeffs(m_coeffs.size());
    std::vector<T> other_coeffs = other.get_coeffs();
    int i_0, j_0, idx_0; //index offsets
    int idx;
    int coeff=0;

    for (int deg0=0; deg0<=m_degree; deg0++){ //loop over order of terms of poly0
        int max_i=m_J[m_nvar][deg0];
        if (deg0==0) i_0=0;
        else i_0=m_N[m_nvar][deg0-1];

        for (int deg1=0; deg1<=m_degree-deg0; deg1++){ //loop over other of terms of poly1
            int max_j=m_J[m_nvar][deg1];
            if (deg1==0) j_0=0;
            else j_0=m_N[m_nvar][deg1-1];

            int deg = deg0+deg1;//order of terms of res
            if (deg==0) idx_0=0;
            else idx_0=m_N[m_nvar][deg-1];

            for (int i=0;i<max_i;i++){
                for (int j=0;j<max_j;j++){
                    if(fabs(m_coeffs[i_0+i])>ZERO && fabs(other_coeffs[j_0+j])>ZERO){
                        if (use_M) idx = m_M[coeff];
                        else{
                        //find what term is the result contributing to
                            std::vector<int> row0=this->get_row(i,deg0);
                            std::vector<int> row1=this->get_row(j,deg1);
                            std::vector<int> row(m_nvar);
                            for (int k=0;k<m_nvar;k++) row[k]=row0[k]+row1[k];
                            idx = res.get_idx(row);
                        }
                        //multiply and add contribution
                        res_coeffs[idx_0+idx]+= m_coeffs[i_0+i]*other_coeffs[j_0+j];
                    }
                    coeff++;
                }
            }

        }
    }

    res.set_coeffs(res_coeffs);
    return res;
}

template <class T>
polynomial<T> polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
polynomial<T> polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
polynomial<T> polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
polynomial<T> polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}


/******************************/
/*UNARY OPERATORS             */
/******************************/
template <class T>
polynomial<T> polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
polynomial<T> polynomial<T>::operator-() const{

    std::vector<T> coeffs=this->get_coeffs();
    polynomial<T> res(m_nvar,m_degree);
    for (int i=0;i<coeffs.size();i++){
        coeffs[i] = -coeffs[i];
    }
    res.set_coeffs(coeffs);
    return res;
}


/******************************/
/*ASSIGNEMENT (with operators)*/
/******************************/
template <class T>
polynomial<T>& polynomial<T>::operator=(const polynomial<T> &other){

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }

    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    m_coeffs = other.get_coeffs();
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator+=(const polynomial<T> &other){
    *this = polynomial<T>::operator+(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator-=(const polynomial<T> &other){
    *this = polynomial<T>::operator-(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator*=(const polynomial<T> &other){
    *this = polynomial<T>::operator*(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator+=(const T& other){
    *this = polynomial<T>::operator+(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator-=(const T& other){
    *this = polynomial<T>::operator-(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator*=(const T& other){
    *this = polynomial<T>::operator*(other);
    return *this;
}

template <class T>
polynomial<T>& polynomial<T>::operator/=(const T& other){
    *this = polynomial<T>::operator/(other);
    return *this;
}

/******************************/
/*COMPARISON                  */
/******************************/
template <class T>
bool polynomial<T>::operator==(const polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(m_coeffs==other.get_coeffs())
        return true;
    return false;
}

template <class T>
bool polynomial<T>::operator!=(const polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(polynomial<T>::operator==(other)) return false;
    else return true;

}

/******************************/
/*EVALUTAION & COMPOSITION    */
/******************************/

//Evaluate canonical base 1, x, x2, x3... in a polynomial.
template <class T>
std::vector<polynomial<T> > polynomial<T>::evaluate_base1D(const polynomial<T> &other) const{

    std::vector<polynomial<T> > v;

    for(int i=0; i<=m_degree; i++){
        v.push_back(polynomial<T>(m_nvar,m_degree));
    }

    v[0] = 1.0;

    for (int i=1; i<=m_degree; i++){
        v[i] = other * v[i-1];
    }
    return v;
}

//1-d Evaluation method
template <class T>
std::vector<T> polynomial<T>::evaluate_base1D(const T &other) const {

    std::vector<T> v(m_degree+1);

    v[0] = 1.0;

    for (int i=1; i<=m_degree; i++){
        v[i] = other * v[i-1];
    }
    return v;
}

//Multivariate Evaluation method
template <class T>
T polynomial<T>::evaluate(const std::vector<T> &x) const {
    if(m_nvar!=x.size()){
        smart_exception(m_name+"(evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }

    //evaluate the bases
    std::vector < std::vector <T> > base;
    base.resize(m_nvar);
    for (int i=0; i<m_nvar;i++){
        base[i].resize(m_degree+1);
        base[i][0]=1.0;
        for (int j=1; j<=m_degree; j++){
            base[i][j]=x[i]*base[i][j-1];
        }
    }

    //construct the full polynomial value
    T res = 0;
    int idx=0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            T prod = 1.0;
            if (fabs(m_coeffs[idx])>ZERO){
                std::vector<int> row = this->get_row(i,deg);
                for(int j=0;j<m_nvar; j++){
                    prod*=base[j][row[j]];
                }
                res += m_coeffs[idx]*prod;
            }
            idx++;
        }
    }

    return res;
}


template <class T>
polynomial<T> polynomial<T>::composition(const std::vector<polynomial<T> > &other) const{
    if(m_nvar!=other.size()){
        smart_exception(m_name+"Composition is with a vector of polynomial of the same size of nvar");
    }

    for(int i=0; i<m_nvar; i++){
        if(m_nvar!=other[i].get_nvar()){
            smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
        }
        if(m_degree!=other[i].get_degree()){
            smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
        }
    }

    //allocate memory
    std::vector<std::vector<polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(polynomial<T>(m_nvar,m_degree));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<polynomial<T> > v = evaluate_base1D(other[j]);
        for(int i=0; i<=m_degree; i++){
            base[j][i] = v[i];
        }
    }

    //composing
    polynomial<T> res(m_nvar,m_degree);
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = this->get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            polynomial<T> prod(m_nvar,m_degree);
            prod.set_coeffs(0,1.0);
            for(int j=0;j<m_nvar; j++){
                prod*=base[j][row[j]];
            }
            res += m_coeffs[count]*prod;
            count++;
        }
    }

    return res;
}



/******************************/
/*J & N matrix manipulation   */
/******************************/
template < class T >
void polynomial<T>::initialize_J(){
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
void polynomial<T>::initialize_N(){
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
std::vector<int> polynomial<T>::get_row(const int &idx, const int &deg) const{

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
int polynomial<T>::get_idx(const std::vector<int> &k) const{
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
void polynomial<T>::initialize_M(const int &nvar, const int &degree){

    polynomial<T> poly(nvar,degree);
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
void polynomial<T>::delete_M(){
    m_M.resize(1);
    m_Mnvar=0;
    m_Mdegree=0;
}


template class polynomial<double>;
template class polynomial<float>;
template class polynomial<long double>;
