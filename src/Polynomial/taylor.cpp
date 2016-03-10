/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "Polynomial/taylor.h"


using namespace smart;
using namespace polynomial;

/******************************/
/*CONSTRUCTORS                */
/******************************/
template < class T >
taylor_polynomial<T>::taylor_polynomial(const int &vars, const int &order) : base_polynomial<T>(vars,order){
    m_name="Taylor Polynomial";
    m_monomial_base = true;
}

template < class T >
taylor_polynomial<T>::taylor_polynomial(const int &vars, const int &order, const int &i) : base_polynomial<T>(vars,order,i){
    m_name="Taylor Polynomial";
    m_monomial_base = true;
}

template < class T >
taylor_polynomial<T>::taylor_polynomial(const int &vars, const int &order, const T &value) : base_polynomial<T>(vars,order,value){
    m_name="Taylor Polynomial";
    m_monomial_base = true;
}

template < class T >
taylor_polynomial<T>::taylor_polynomial(const int &vars, const int &order, const int &i, const T &a, const T &b) : base_polynomial<T>(vars,order,i,a,b){
    m_name="Taylor Polynomial";
    m_monomial_base = true;
}


template < class T >
taylor_polynomial<T>::~taylor_polynomial(){
}

/******************************/
/*ARITHMETIC OPERATIONS (+-*) */
/******************************/
template < class T >
taylor_polynomial<T> taylor_polynomial<T>::operator+(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    taylor_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;


}

template < class T >
taylor_polynomial<T> taylor_polynomial<T>::operator-(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    taylor_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}


template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator*(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    taylor_polynomial<T> res(m_nvar,m_degree);
    base_polynomial<T>::monomial_multiplication(*this,other,res);

    return res;
}

template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator/(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    taylor_polynomial<T> res(m_nvar,m_degree);
    res = inv(other);

    return res*(*this);
}


template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}


/******************************/
/*UNARY OPERATORS             */
/******************************/
template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
taylor_polynomial<T> taylor_polynomial<T>::operator-() const{

    std::vector<T> coeffs=this->get_coeffs();
    taylor_polynomial<T> res(m_nvar,m_degree);
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
taylor_polynomial<T>& taylor_polynomial<T>::operator=(const taylor_polynomial<T> &other){

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    m_coeffs = other.get_coeffs();
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator+=(const taylor_polynomial<T> &other){
    *this = taylor_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator-=(const taylor_polynomial<T> &other){
    *this = taylor_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator*=(const taylor_polynomial<T> &other){
    *this = taylor_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator/=(const taylor_polynomial<T> &other){
    *this = taylor_polynomial<T>::operator/(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator+=(const T& other){
    *this = taylor_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator-=(const T& other){
    *this = taylor_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator*=(const T& other){
    *this = taylor_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
taylor_polynomial<T>& taylor_polynomial<T>::operator/=(const T& other){
    *this = taylor_polynomial<T>::operator/(other);
    return *this;
}


template < class T >
taylor_polynomial<T> taylor_polynomial<T>::inv(const taylor_polynomial<T> &other) const{

         std::vector <T> coeffs = this -> get_coeffs(); // p = c + n(x)
         T c = coeffs[0];

         if (fabs(c)<ZERO){
             smart_throw(m_name+": Division by zero occurred");
         }

         taylor_polynomial<T> times(m_nvar,m_degree); // times = -n/c
         coeffs[0] = 0.0;
         times.set_coeffs(coeffs);
         times/= -c;

         taylor_polynomial<T> res(m_nvar,m_degree, (T) (1/c));
         taylor_polynomial<T> term = res;

         for (int i=1; i<=m_degree; i++){ // 1/p = 1/c - n/c^2 + n^2/c^3 + ...
             term*=times;
             res+=term;
         }

         return res;
}


/******************************/
/*COMPARISON                  */
/******************************/
template <class T>
bool taylor_polynomial<T>::operator==(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(m_coeffs==other.get_coeffs())
        return true;
    return false;
}

template <class T>
bool taylor_polynomial<T>::operator!=(const taylor_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(taylor_polynomial<T>::operator==(other)) return false;
    else return true;

}


/******************************/
/*EVALUATION & COMPOSITION    */
/******************************/

//evaluate chebyshev base t0(x), t1(x), t2(x) in a polynomial. It first map x from [a,b] to [-1,1]
template <class T>
std::vector<taylor_polynomial<T> > taylor_polynomial<T>::evaluate_base1D(const taylor_polynomial<T> &other){
    int nvar = other.get_nvar();
    int degree = other.get_degree();

    std::vector<taylor_polynomial<T> > v;

    for(int i=0; i<=degree; i++){
        v.push_back(taylor_polynomial<T>(nvar,degree));
    }

    for(int i=0;i<=degree;i++)
        base_polynomial<T>::evaluate_base1D_monomial(i,other,v[i]);

    return v;
}

template <class T>
void taylor_polynomial<T>::composition(const std::vector<taylor_polynomial<T> > &other) {
    if(m_nvar!=other.size()){
        smart_throw(m_name+": Composition is with a vector of polynomial of the same size of nvar");
    }

    for(int i=0; i<m_nvar; i++){
        if(m_monomial_base != other[i].is_monomial_base()){
            smart_throw(m_name+": One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
        }
        if(m_nvar!=other[i].get_nvar()){
            smart_throw(m_name+": Polynomials don't have the same number of variables. They don't belong to the same Algebra");
        }
        if(m_degree!=other[i].get_degree()){
            smart_throw(m_name+": Polynomials don't have the same order. They don't belong to the same Algebra");
        }
    }

    //allocate memory
    std::vector<std::vector<taylor_polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<taylor_polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(taylor_polynomial<T>(m_nvar,m_degree));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<taylor_polynomial<T> > v = evaluate_base1D(other[j]);
        for(int i=0; i<=m_degree; i++){
            base[j][i] = v[i];
        }
    }

    //composing
    taylor_polynomial<T> res(m_nvar,m_degree);
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = this->get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            taylor_polynomial<T> prod(m_nvar,m_degree);
            prod.set_coeffs(0,1.0);
            for(int j=0;j<m_nvar; j++){
                prod*=base[j][row[j]];
            }
            res += m_coeffs[count]*prod;
            count++;
        }
    }

    (*this) = res;
}

template <class T>
std::vector<T> taylor_polynomial<T>::evaluate_basis(const std::vector<T> &x) const { //most direct implementation, faster ones might be available
    return this->evaluate_basis_monomial(x);
}

//Multivariate Evaluation method
template <class T>
T taylor_polynomial<T>::evaluate(const std::vector<T> &x) const { //most direct implementation, faster ones might be available

    std::vector<T> basis = evaluate_basis(x);

    //construct the full polynomial value
    T res = 0;
    for(int i=0;i<m_coeffs.size(); i++)
        res+=m_coeffs[i]*basis[i];

    return res;

}

//1-d Evaluation method
template <class T>
T taylor_polynomial<T>::evaluate(const T &x) const {
    if(m_nvar>1){
        smart_throw(m_name+": (evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }

    return m_coeffs[0]+this->horner(x,1);

}


template < class T >
void taylor_polynomial<T>::map(const int &idx, const std::vector<T> &a, const std::vector<T> &b){

    if(b.size() != a.size())
        smart_throw(m_name+": mapping of polynomial variable from [-1,1]^d to [a,b]^d a and b need to be vector of the same size");

    std::vector<taylor_polynomial<T> > mapped_vars;

    // construct polynomial, x1, x2, x3,...
    for(int i=0; i<m_nvar; i++){
        if(b[i]<=a[i])
            smart_throw(m_name+": mapping of polynomial variable from [-1,1] to [a,b] with b<=a");
        mapped_vars.push_back(taylor_polynomial<T>(m_nvar, m_degree,i));
        mapped_vars[i] = (b[i]-a[i])/2.0 * mapped_vars[i] + (b[i]+a[i])/2.0;
    }

    composition(mapped_vars);

}


/******************************/
/*BASIS MANIPULATION          */
/******************************/

template < class T >
void taylor_polynomial<T>::to_monomial_basis(){
    if(m_monomial_base)
        smart_throw(m_name+": The transformation to monomial bases has been called when the base is already monomial.");
}


template < class T >
void taylor_polynomial<T>::from_monomial_basis(){
    if(m_monomial_base)
        smart_throw(m_name+": The transformation from monomial bases has been called when the base is not monomial.");
}

template < class T >
std::string taylor_polynomial<T>::get_basis_name() const{
    return "T";
}

/******************************/
/*APPROXIMATION               */
/******************************/
template < class T >
std::vector<T> taylor_polynomial<T>::approximation(T (*f)(T x), const T &x0){

    return std::vector<T>();
}



template class taylor_polynomial<double>;
template class taylor_polynomial<float>;
template class taylor_polynomial<long double>;

