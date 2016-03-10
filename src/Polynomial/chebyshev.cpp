/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "Polynomial/chebyshev.h"

using namespace smart;
using namespace polynomial;

/******************************/
/*CONSTRUCTORS                */
/******************************/
template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const bool& monomial) : base_polynomial<T>(vars,order){
    m_name="Chebyshev Polynomial";
    m_monomial_base=monomial;
}

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const int &i, const bool& monomial) : base_polynomial<T>(vars,order,i){
    m_name="Chebyshev Polynomial";
    m_monomial_base=monomial;
}

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const T &value, const bool& monomial) : base_polynomial<T>(vars,order,value){
    m_name="Chebyshev Polynomial";
    m_monomial_base=monomial;
}

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const int &i, const T &a, const T &b, const bool& monomial) : base_polynomial<T>(vars,order,i,a,b){
    m_name="Chebyshev Polynomial";
    m_monomial_base=monomial;
}

template < class T >
chebyshev_polynomial<T>::~chebyshev_polynomial(){

}

/******************************/
/*ARITHMETIC OPERATIONS (+-*) */
/******************************/
template < class T >
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    chebyshev_polynomial<T> res(m_nvar,m_degree,other.is_monomial_base());

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;


}

template < class T >
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    chebyshev_polynomial<T> res(m_nvar,m_degree,other.is_monomial_base());

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}


//OPERATOR* OVERLOADING FOR DCT-BASED MULTIPLICATION
//Author : Carlos Ortega Absil (carlos.ortega@strath.ac.uk)
//Note: the indexing, scaling, etc. operations could be suppressed with 1-var polynomials for a performance gain
template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator*(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    //perform multiplication in monomial base
    if(m_monomial_base){
        chebyshev_polynomial<T> res(m_nvar,m_degree,other.is_monomial_base());
        base_polynomial<T>::monomial_multiplication(*this,other,res);
        return res;
    }
    else
    {
    #ifdef CHEBYSHEV_DCT_MULTIPLICATION
    int ncoeffs = combination(m_nvar,m_degree);
    std::vector<T> other_coeffs = other.get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    //initialise stuff needed by fftw
    int dct_degree_aux = int (m_degree*1.5+1);
    int dct_degree[m_nvar];
    for (int i=0;i<m_nvar;i++){
        dct_degree[i]=dct_degree_aux+1;
    }
    int pointers_length=pow(dct_degree_aux+1,m_nvar);
    T *dct0, *dct1, *dct01;

    // allocate pointers for all DCTs
    dct_malloc(dct0,pointers_length);
    dct_malloc(dct1,pointers_length);
    dct_malloc(dct01,pointers_length);

    // build a vector idx[ncoeffs] with the pointer indexes in row-major format
    // build a vector scale[ncoeffs] with the scale factor for each of the coefficients
    std::vector <int> idx;
    std::vector <T> scale;
    int ii=0;
    for (int deg=0;deg<=m_degree;deg++){
        int max_i=m_J[m_nvar][deg];
        for (int i=0;i<max_i;i++){
            std::vector<int> row = res.get_row(i,deg);
            int term_idx=0;
            T term_scale=1.0;
            for (int var=0;var<m_nvar;var++){
                if (row[var]==0) term_scale*=2.0;
                term_idx+=row[var]*pow((dct_degree_aux+1),var);
            }
            idx.push_back(term_idx);
            scale.push_back(term_scale);
            // initialise non-zero terms of dct0 and dct1 here too (to avoid an additional for loop)
            dct0[term_idx]=m_coeffs[ii]*term_scale;
            dct1[term_idx]=other_coeffs[ii]*term_scale;
            ii++;
        }
    }

    //DCT(x0) and DCT(x1)
    dct_do(m_nvar,dct_degree,dct0);
    dct_do(m_nvar,dct_degree,dct1);

    // component-wise multiplication DCT(x0):DCT(x1)
    // scale already to avoid coefficients growing too much in large algebras
    T scale_intermediate = pow(4*dct_degree_aux,m_nvar);
    for(int i=0;i<pointers_length;i++){
        dct01[i]=dct0[i]*dct1[i]/scale_intermediate;
    }

    // deallocate more stuff
    dct_free(dct0);
    dct_free(dct1);

    // Obtain DCT(DCT(x0):DCT(x1))
    dct_do(m_nvar,dct_degree,dct01);

    // rescale and save results
    for (int i=0;i<ncoeffs;i++){
        T term_result=dct01[idx[i]]/scale[i];
        if (fabs(term_result)>ZERO) res.set_coeffs(i,term_result);
        else res.set_coeffs(i,0.0);
    }

    // deallocate and return
    dct_free(dct01);
    return res;

    #else

    return direct_multiplication(*this,other);

    #endif
    }
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator/(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    chebyshev_polynomial<T> res(m_nvar,m_degree,other.is_monomial_base());
    res = inv(other);

    return res*(*this);
}


template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}


/******************************/
/*UNARY OPERATORS             */
/******************************/
template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-() const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree,m_monomial_base);
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
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator=(const chebyshev_polynomial<T> &other){

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    m_coeffs = other.get_coeffs();
    m_monomial_base=other.is_monomial_base();
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator+=(const chebyshev_polynomial<T> &other){
    *this = chebyshev_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator-=(const chebyshev_polynomial<T> &other){
    *this = chebyshev_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator*=(const chebyshev_polynomial<T> &other){
    *this = chebyshev_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator/=(const chebyshev_polynomial<T> &other){
    *this = chebyshev_polynomial<T>::operator/(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator+=(const T& other){
    *this = chebyshev_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator-=(const T& other){
    *this = chebyshev_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator*=(const T& other){
    *this = chebyshev_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator/=(const T& other){
    *this = chebyshev_polynomial<T>::operator/(other);
    return *this;
}


template < class T >
chebyshev_polynomial<T> chebyshev_polynomial<T>::inv(const chebyshev_polynomial<T> &other) const{
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree,other.is_monomial_base());

    std::vector<T> range = other.get_range();
    T a, b;

    a = range[0];
    b = range[1];
    std::vector<T> cheb_inv = approximation(inverse,a,b);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = evaluate_base1D(other, a,b);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_inv[i];
    }

    return res;
}

// DIRECT MULTIPLICATION
template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T> :: direct_multiplication(const chebyshev_polynomial<T> &x0, const chebyshev_polynomial<T> &x1) const{
    if(x0.get_nvar()!=x1.get_nvar()){
        smart_exception("Direct multiplication: Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(x0.get_degree()!=x1.get_degree()){
        smart_exception("Direct multiplication: Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    chebyshev_polynomial<T> res(x0.get_nvar(),x0.get_degree());
    std::vector<T> res_coeffs(combination(x0.get_nvar(),x0.get_degree()));
    double nvariations = pow(2,x0.get_nvar());
    std::vector<T> x0_coeffs = x0.get_coeffs();
    std::vector<T> x1_coeffs = x1.get_coeffs();
    std::vector<std::vector<int> > x0_J=x0.get_J();
    std::vector<std::vector<int> > x0_N=x0.get_N();
    res.initialize_t();
    std::vector<std::vector<int> > t=res.get_t();
    for(int i=0; i<=x0.get_degree(); i++){//loop over subset degree i of poly1
        for(int j=0; j<=x1.get_degree(); j++){//loop over subset degree j of poly2
            //if((i+j)<=m_degree){
                for(int idx1=0; idx1<x0_J[x0.get_nvar()][i]; idx1++){//index over elements with degree i in poly1
                    for(int idx2=0; idx2<x0_J[x0.get_nvar()][j]; idx2++){//index over elements with degree j in poly2
                        int sub_idx1=0, sub_idx2=0, sub_idx3=0;
                        if(i>0) sub_idx2=x0_N[x0.get_nvar()][i-1];
                        if(j>0) sub_idx3=x0_N[x0.get_nvar()][j-1];
                        if(fabs(x0_coeffs[sub_idx2+idx1])>ZERO && fabs(x1_coeffs[sub_idx3+idx2])>ZERO){
                            std::vector<int> v1 = x0.get_row(idx1,i);
                            std::vector<int> v2 = x0.get_row(idx2,j);
                            std::vector<int> v3(x0.get_nvar());
                            T term = (1.0/nvariations)*(x0_coeffs[sub_idx2+idx1]*x1_coeffs[sub_idx3+idx2]);
                            for(int iter=0; iter<nvariations; iter++){
                                for(int k=0; k<x0.get_nvar(); k++){
                                    v3[k] = std::fabs(v1[k]+t[iter][k]*v2[k]);
                                }
                                int deg3 = std::accumulate(v3.begin(),v3.end(),0);
                                if(deg3<=x0.get_degree()){
                                    int pos = res.get_idx(v3);
                                    sub_idx1 = 0;
                                    if(deg3>0) sub_idx1=x0_N[x0.get_nvar()][deg3-1];
                                    res_coeffs[sub_idx1 + pos] += term;
                                }
                            }
                        }
                    }
                }
            //}
        }
    }

    res.set_coeffs(res_coeffs);
    return res;
}



/******************************/
/*COMPARISON                  */
/******************************/
template <class T>
bool chebyshev_polynomial<T>::operator==(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(m_coeffs==other.get_coeffs() && m_monomial_base==other.is_monomial_base())
        return true;
    return false;
}

template <class T>
bool chebyshev_polynomial<T>::operator!=(const chebyshev_polynomial<T> &other) const{

    if(m_monomial_base != other.is_monomial_base()){
        smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
    }
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(chebyshev_polynomial<T>::operator==(other)) return false;
    else return true;

}


/******************************/
/*EVALUATION & COMPOSITION    */
/******************************/

//evaluate chebyshev base t0(x), t1(x), t2(x) in a polynomial. It first map x from [a,b] to [-1,1]
template <class T>
std::vector<chebyshev_polynomial<T> > chebyshev_polynomial<T>::evaluate_base1D(const chebyshev_polynomial<T> &other, const T &a, const T &b) {
    if(b<a)
        smart_exception("Base evaluation is in range [a,b] where b<a");

    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector<chebyshev_polynomial<T> > v;

    for(int i=0; i<=degree; i++){
        v.push_back(chebyshev_polynomial<T>(nvar,degree,other.is_monomial_base()));
    }
    
    chebyshev_polynomial<T> mapped(nvar,degree,other.is_monomial_base());
    if(b==a)
        mapped.set_coeffs(0,a);
    else
        mapped = (2.0*(other)-(a+b))/(b-a);

    v[0] = 1.0;

    if(other.is_monomial_base()){
        for (int i=1; i<=degree; i++){
            v[i] = mapped * v[i-1];
        }
        return v;
    }

    else{
        v[1] = mapped;

        for (int i=2; i<=degree; i++){
            v[i] = 2.0 * mapped * v[i-1] - v[i-2];
        }

    }
    return v;
}

template <class T>
void chebyshev_polynomial<T>::composition(const std::vector<chebyshev_polynomial<T> > &other){
    if(m_nvar!=other.size()){
        smart_exception(m_name+"Composition is with a vector of polynomial of the same size of nvar");
    }

    for(int i=0; i<m_nvar; i++){
        if(m_monomial_base != other[i].is_monomial_base()){
            smart_exception(m_name+"One of the two polynomials has not been transformed to monomial base. They do not belong to the same Algebra");
        }
        if(m_nvar!=other[i].get_nvar()){
            smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
        }
        if(m_degree!=other[i].get_degree()){
            smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
        }
    }

    //allocate memory
    std::vector<std::vector<chebyshev_polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<chebyshev_polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(chebyshev_polynomial<T>(m_nvar,m_degree,other[j].is_monomial_base()));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<chebyshev_polynomial<T> > v = evaluate_base1D(other[j],-1.0,1.0);
        for(int i=0; i<=m_degree; i++)
            base[j][i] = v[i];
    }

    //composing
    chebyshev_polynomial<T> res(m_nvar,m_degree,other[0].is_monomial_base());
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = this->get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            chebyshev_polynomial<T> prod(m_nvar,m_degree,other[0].is_monomial_base());
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
std::vector<T> chebyshev_polynomial<T>::evaluate_basis(const std::vector<T> &x) const{
    if(m_nvar!=x.size()){
        smart_exception(m_name+"(evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }
    for (int i=0;i<m_nvar;i++){
        if (fabs(x[i])>1){
            smart_exception(m_name+"(evaluate) All components of point must belong to [-1,1].");
        }
    }


    if(m_monomial_base){
        return this->evaluate_basis_monomial(x);
    }

    //evaluate the bases
    std::vector < std::vector <T> > base;
    base.resize(m_nvar);
    for (int i=0; i<m_nvar;i++){
        base[i].resize(m_degree+1);
        base[i][0]=1.0;
        if (m_degree>0) base[i][1]=x[i];
        for (int j=2; j<=m_degree; j++){
            base[i][j]=2.0*x[i]*base[i][j-1]-base[i][j-2];
        }
    }


    //evaluate the full polynomial bases
    std::vector<T> res(m_coeffs.size());
    int idx=0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            T prod = 1.0;
            if (fabs(m_coeffs[idx])>ZERO){
                std::vector<int> row = this->get_row(i,deg);
                for(int j=0;j<m_nvar; j++){
                    prod*=base[j][row[j]];
                }
                res[idx] = prod;
            }
            idx++;
        }
    }

    return res;
}

//Multivariate Evaluation method
//Author: Carlos Ortega Absil (carlos.ortega@strath.ac.uk)
template <class T>
T chebyshev_polynomial<T>::evaluate(const std::vector<T> &x) const { //most direct implementation, faster ones might be available

    std::vector<T> basis = evaluate_basis(x);

    //construct the full polynomial value
    T res = 0;
    for(int i=0;i<m_coeffs.size(); i++)
        res+=m_coeffs[i]*basis[i];

    return res;
}

//1-d Evaluation method
//Author: Carlos Ortega Absil (carlos.ortega@strath.ac.uk)
template <class T>
T chebyshev_polynomial<T>::evaluate(const T &x) const {
    if(m_nvar>1){
        smart_exception(m_name+"(evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }
    if (fabs(x)>1){
        smart_exception(m_name+"(evaluate) All components of point must belong to [-1,1].");
    }

    if(m_monomial_base){
        return m_coeffs[0]+this->horner(x,1);
    }

    return m_coeffs[0]+x*clenshaw(x,1)-clenshaw(x,2);
}

/******************************/
/*BASIS MANIPULATION          */
/******************************/

template < class T >
void chebyshev_polynomial<T>::to_monomial_basis(){

    if(m_monomial_base)
        smart_exception(m_name+"The transformation to monomial bases has been called when the base is already monomial.");

    m_monomial_base=true;

    //check that polynomial is not constant neither 1st degree. In this case do nohing
    T sum = 0.0;
    for(int i=m_nvar+1; i<m_coeffs.size(); i++)
        sum += fabs(m_coeffs[i]);
    if(sum==0)
        return;

    chebyshev_polynomial<T> res(m_nvar,m_degree,(T) 0.0, true);

    int ncoeffs=res.get_coeffs().size();

    std::vector <chebyshev_polynomial <T> > term_vector;

    for (int i=0;i<ncoeffs;i++){
        chebyshev_polynomial<T> coeff(m_nvar,m_degree,(T) m_coeffs[i],true);
        term_vector.push_back(coeff);
    }

    for (int v=0;v<m_nvar;v++){
        chebyshev_polynomial<T> base2(m_nvar,m_degree,(T) 1.0, true);
        chebyshev_polynomial<T> base1(m_nvar,m_degree,(int) v, true);
        chebyshev_polynomial<T> x(m_nvar,m_degree,(int) v, true);
        chebyshev_polynomial<T> term(m_nvar,m_degree, true);

        for (int d=1;d<=m_degree;d++){
            if (d==1)  term=base1;
            else{
                term=2.0*x*base1-base2;
                base2=base1;
                base1=term;
            }
            //term is the chebyshev term of order d in variable v. Now we multiply by it the necessary terms in term_vector
            int coeff_idx=m_N[m_nvar][d-1];
            for(int deg=d; deg<=m_degree; deg++){
                for(int i=0; i<m_J[m_nvar][deg]; i++){
                    std::vector<int> row = res.get_row(i,deg);
                    if (row[v]==d) term_vector[coeff_idx]*=term;
                    coeff_idx+=1;
                }
            }
        }
    }

    for (int i=0;i<ncoeffs;i++){
        res+=term_vector[i];
    }

    m_coeffs = res.get_coeffs();

}


template < class T >
void chebyshev_polynomial<T>::from_monomial_basis(){
    if(!m_monomial_base)
        smart_exception(m_name+"The transformation from monomial bases has been called when the base is not in monomial.");

    m_monomial_base=false;

    chebyshev_polynomial<T> res(m_nvar,m_degree,(T) 0.0);
    int ncoeffs=res.get_coeffs().size();

    std::vector <chebyshev_polynomial <T> > term_vector;

    for (int i=0;i<ncoeffs;i++){
        term_vector.push_back(chebyshev_polynomial<T>(m_nvar,m_degree,(T) m_coeffs[i]));
    }

    for (int v=0;v<m_nvar;v++){
        chebyshev_polynomial<T> base(m_nvar,m_degree, (T) 1.0);
        chebyshev_polynomial<T> x(m_nvar,m_degree,(int) v);
        chebyshev_polynomial<T> term(m_nvar,m_degree);
        for (int d=1;d<=m_degree;d++){
            term=x*base;
            base=term;
            //term is the monomial of order d in variable v (aka v^d), in chebyshev basis. Now we multiply by it the necessary terms in term_vector
            int coeff_idx=m_N[m_nvar][d-1];
            for(int deg=d; deg<=m_degree; deg++){
                for(int i=0; i<m_J[m_nvar][deg]; i++){
                    std::vector<int> row = res.get_row(i,deg);
                    if (row[v]==d) term_vector[coeff_idx]*=term;
                    coeff_idx+=1;
                }
            }
        }
    }

    for (int i=0;i<ncoeffs;i++){
        res+=term_vector[i];
    }

    m_coeffs = res.get_coeffs();

}

template < class T >
std::string chebyshev_polynomial<T>::get_basis_name() const{
    return "C";
}

template < class T >
void chebyshev_polynomial<T>::map(const int &idx, const std::vector<T> &a, const std::vector<T> &b){

    if(b.size() != a.size())
        smart_exception(m_name+"mapping of polynomial variable from [-1,1]^d to [a,b]^d a and b need to be vector of the same size");

    std::vector<chebyshev_polynomial<T> > mapped_vars;

    // construct polynomial, x1, x2, x3,...
    for(int i=0; i<m_nvar; i++){
        if(b[i]<=a[i])
            smart_exception(m_name+"mapping of polynomial variable from [-1,1] to [a,b] with b>=a");
        mapped_vars.push_back(chebyshev_polynomial<T>(m_nvar, m_degree,i, m_monomial_base));
        mapped_vars[i] = (b[i]-a[i])/2.0 * mapped_vars[i] + (b[i]+a[i])/2.0;
    }

    composition(mapped_vars);

}

/******************************/
/*APPROXIMATION               */
/******************************/
template < class T >
std::vector<T> chebyshev_polynomial<T>::approximation(T (*f)(T x), const T &a, const T &b, const T &deg){
    // returns a vector of size deg+1 with the coefficients of the chebyshev approximation of an univariate function
    int n = chebyshev_polynomial<T>::MAX_DEGREE;
    std::vector<T> res(deg+1), d(n+1);
    T fac;
    T pi = 3.141592653589793;
    T t;
    T total;
    T y;

    for (int k = 0; k <= n; k++)
    {
        t = cos(pi*(k+0.5)/(n+1)); //zeros Ti
        y = ((1.0+t)*b + (1.0-t)*a)/2.0; //mapped zeros
        d[k] = f(y); //evaluate function
    }

    //Interpolation in chebyshev basis -> MAX_DEGREE chebyshev nodes but only a deg polynomial.
    fac = 2.0/(n+1);
    for (int j = 0; j <= deg; j++)
    {
        total = 0.0;
        for (int k = 0; k <= n; k++)
        {
            total = total+d[k]*cos( (pi*j)*( (k+ 0.5)/(n+1) ) );
        }
        res[j] = fac*total;
    }

    res[0] = res[0]/2.0;

    return res;

}

template < class T >
chebyshev_polynomial<T> chebyshev_polynomial<T>::approximation(T (*f)(T x), const chebyshev_polynomial<T> &other){
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree, other.is_monomial_base());
    std::vector<T> range = other.get_range();
    std::vector<T> approx(degree+1);


   if (other.is_monomial_base()){

       //approximate sin in [a,b] with increased degree (two-step truncation to enhance precision)
       int deg_max = chebyshev_polynomial<T>::MAX_DEGREE;
       int deg = std::min((int) (degree*1.5+1), deg_max);
       std::vector<T> cheb_approx = chebyshev_polynomial<T>::approximation(f,range[0],range[1],deg);

       // Translation to canonical basis, taking into acount deg+1 terms from cheb_approx but building a monom_approx of degree+1 terms.
       // Hence rewriting code instead of calling to_monomial(), to avoid the computation of worthless terms of order > degree.

       chebyshev_polynomial<T> monom_approx(1,degree,(T) cheb_approx[0], true);
       chebyshev_polynomial<T> x(1,degree,(int) 0, true);
       chebyshev_polynomial<T> cheb_base1(1,degree,(int) 0, true);
       chebyshev_polynomial<T> cheb_base2(1,degree, (T) 1.0, true);
       chebyshev_polynomial<T> cheb_base(1,degree, true);

       monom_approx+= cheb_approx[1]*x;
       for (int i=2;i<=deg;i++){
           cheb_base=2.0*x*cheb_base1-cheb_base2;
           cheb_base2=cheb_base1;
           cheb_base1=cheb_base;
           monom_approx+=cheb_approx[i]*cheb_base;
       }

       approx = monom_approx.get_coeffs();
   }
   else{
       //approximate sin in [a,b]
       approx = chebyshev_polynomial<T>::approximation(f,range[0],range[1],degree);
   }

    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = evaluate_base1D(other,range[0],range[1]);
    for (int i=0; i<=degree; i++){
        res += base[i]*approx[i];
    }

    return res;
}

/******************************/
/*PRIVATE ROUTINES            */
/******************************/
template < class T >
void chebyshev_polynomial<T>::initialize_t(){

    m_t.resize(pow(2,m_nvar));
    for(int i=0; i<m_nvar; i++){
        m_t[i].resize(2);
    }

    std::vector<int> values(2);
    values[0] = -1;
    values[1] = 1;

    variations(values,m_nvar, m_t);

}

// private routine for evaluation
template <class T>
T chebyshev_polynomial<T>::clenshaw(T x, int n) const{
    if (n>m_degree) return 0;
    else return m_coeffs[n]+2*x*clenshaw(x,n+1)-clenshaw(x,n+2);
}
 
template class chebyshev_polynomial<double>;
template class chebyshev_polynomial<float>;
template class chebyshev_polynomial<long double>;
