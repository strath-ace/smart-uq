/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/


#include "Polynomial/chebyshev.h"


using namespace smart;
using namespace polynomial;

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order) : base_polynomial<T>("Chebyshev Polynomial", vars,order){

}

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const int &i)  : base_polynomial<T>("Chebyshev Polynomial", vars,order,i){

}

template < class T >
chebyshev_polynomial<T>::chebyshev_polynomial(const int &vars, const int &order, const T &value) : base_polynomial<T>("Chebyshev Polynomial", vars,order,value){

}


template < class T >
std::string chebyshev_polynomial<T>::get_basis_name() const
{
    return "T";
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+(const chebyshev_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }

    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-(const chebyshev_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

// //OPERATOR* OVERLOADING FOR DIRECT MULTIPLICATION
// template <class T>
// chebyshev_polynomial<T> chebyshev_polynomial<T>::operator*(const chebyshev_polynomial<T> &other) const{
    
//     if(m_nvar!=other.get_nvar()){
//         std::cout<<"Polynomials don't have the same number of variables. They don't belong to the same Algebra"<<std::endl;
//         exit(EXIT_FAILURE);
//     }
//     if(m_degree!=other.get_degree()){
//         std::cout<<"Polynomials don't have the same order. They don't belong to the same Algebra"<<std::endl;
//         exit(EXIT_FAILURE);
//     }

//     chebyshev_polynomial<T> res(m_nvar,m_degree);
//     std::vector<T> res_coeffs(combination(m_nvar,m_degree));
//     double nvariations = pow(2,m_nvar);
//     std::vector<T> other_coeffs = other.get_coeffs();
//     for(int i=0; i<=m_degree; i++){//loop over subset degree i of poly1
//         for(int j=0; j<=other.get_degree(); j++){//loop over subset degree j of poly2
//             //if((i+j)<=m_degree){
//                 for(int idx1=0; idx1<m_J[m_nvar][i]; idx1++){//index over elements with degree i in poly1
//                     for(int idx2=0; idx2<m_J[m_nvar][j]; idx2++){//index over elements with degree j in poly2
//                         int sub_idx1=0, sub_idx2=0, sub_idx3=0;
//                         if(i>0) sub_idx2=m_N[m_nvar][i-1];
//                         if(j>0) sub_idx3=m_N[m_nvar][j-1];
//                         if(fabs(m_coeffs[sub_idx2+idx1])>ZERO && fabs(other_coeffs[sub_idx3+idx2])>ZERO){
//                             std::vector<int> v1 = this->get_row(idx1,i);
//                             std::vector<int> v2 = this->get_row(idx2,j);
//                             std::vector<int> v3(m_nvar);
//                             for(int iter=0; iter<nvariations; iter++){
//                                 for(int k=0; k<m_nvar; k++){
//                                     v3[k] = std::fabs(v1[k]+m_t[iter][k]*v2[k]);
//                                 }
//                                 int deg3 = std::accumulate(v3.begin(),v3.end(),0);
//                                 if(deg3<=m_degree){
//                                     int pos = res.get_idx(v3);
//                                     sub_idx1 = 0;
//                                     if(deg3>0) sub_idx1=m_N[m_nvar][deg3-1];
//                                     res_coeffs[sub_idx1 + pos] +=
//                                         (1.0/nvariations)*(m_coeffs[sub_idx2+idx1]*other_coeffs[sub_idx3+idx2]);
//                                 }
//                             }
//                         }
//                     }
//                 }
//             //}
//         }
//     }

//     res.set_coeffs(res_coeffs);
//     return res;
// }

//OPERATOR* OVERLOADING FOR DCT-BASED MULTIPLICATION
//Author : Carlos Ortega Absil (carlos.ortega@strath.ac.uk)
//Note: the indexing, scaling, etc. operations could be suppressed with 1-var polynomials for a performance gain
template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator*(const chebyshev_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

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


template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::inv(const chebyshev_polynomial<T> &other) const{
    int nvar =  other.get_nvar();
    int degree = other.get_degree();
    chebyshev_polynomial<T> res(nvar,degree);

    //chebyshev expansion of sin in [a,b]
    //    std::vector<T> coeffs = other.get_coeffs();
    std::vector<T> range = other.get_range();
    T a, b;

    //    //computing value of other in zero
    //    T value = 0.0;
    //    int count = 0;
    //    for(int deg=0; deg<=degree; deg++){
    //        for(int i=0; i<m_J[nvar][deg]; i++){
    //            std::vector<int> row = get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
    //            int prod = 1;
    //            for(int j=0;j<m_nvar; j++){
    //                if(row[j]%2 == 0){
    //                    if(row[j]%4 != 0)
    //                        prod*= -1.0;
    //                }
    //                else prod *= 0.0;
    //            }
    //            value += coeffs[count]*prod;
    //            count++;
    //        }
    //    }

    //    T one = 1.0;
    //    if(value>0){
    //        a = min(one,range[1]);
    //        b = max(one,range[1]);
    //    }
    //    else if(value<0){
    //        a = max(-one,range[0]);
    //        b = min(-one,range[0]);
    //    }
    //    else{
    //        std::cout<<"Inverting a function with a zero in the interval [-1,1]."<<std::endl;
    //        exit(EXIT_FAILURE);
    //    }

    a = range[0];
    b = range[1];
    std::vector<T> cheb_inv = cheb_approximation(inverse,a,b);
    //univariate composition
    std::vector<chebyshev_polynomial<T> > base = evaluate_base(other, a,b);
    for (int i=0; i<=degree; i++){
        res += base[i]*cheb_inv[i];
    }

    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator/(const chebyshev_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    chebyshev_polynomial<T> res(m_nvar,m_degree);
    res = inv(other);

    return res*(*this);
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::operator-() const{
    
    std::vector<T> coeffs=this->get_coeffs();
    chebyshev_polynomial<T> res(m_nvar,m_degree);
    for (int i=0;i<coeffs.size();i++){
        coeffs[i]= -coeffs[i];
    }
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
chebyshev_polynomial<T>& chebyshev_polynomial<T>::operator=(const chebyshev_polynomial<T> &other){

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

template <class T>
bool chebyshev_polynomial<T>::operator==(const chebyshev_polynomial<T> &other) const{

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
bool chebyshev_polynomial<T>::operator!=(const chebyshev_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(chebyshev_polynomial<T>::operator==(other)) return false;
    else return true;

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

    return m_coeffs[0]+x*clenshaw(x,1)-clenshaw(x,2);
}

//Multivariate Evaluation method
//Author: Carlos Ortega Absil (carlos.ortega@strath.ac.uk)
template <class T>
T chebyshev_polynomial<T>::evaluate(const std::vector<T> &x) const { //most direct implementation, faster ones might be available
    if(m_nvar!=x.size()){
        smart_exception(m_name+"(evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }
    for (int i=0;i<m_nvar;i++){
        if (fabs(x[i])>1){
            smart_exception(m_name+"(evaluate) All components of point must belong to [-1,1].");
        }
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

//evaluate chebyshev base t0(x), t1(x), t2(x) in a polynomial. It first map x from [a,b] to [-1,1]
template <class T>
std::vector<chebyshev_polynomial<T> > chebyshev_polynomial<T>::evaluate_base(const chebyshev_polynomial<T> &other, const T &a, const T &b){

    int nvar =  other.get_nvar();
    int degree = other.get_degree();

    std::vector<chebyshev_polynomial<T> > v;

    for(int i=0; i<=degree; i++){
        v.push_back(chebyshev_polynomial<T>(nvar,degree));
    }

    //mapping the argument from the domain of composition to [-1,1]
    chebyshev_polynomial<T> mapped(nvar,degree);
    if(b==a)
        mapped.set_coeffs(0,a);
    else
        mapped = (2.0*other-(a+b))/(b-a);

    v[0] = 1.0;
    v[1] = mapped;

    for (int i=2; i<=degree; i++){
        v[i] = 2.0 * mapped * v[i-1] - v[i-2];
    }

    return v;
}

template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T>::composition(const std::vector<chebyshev_polynomial<T> > &other) const{
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
    std::vector<std::vector<chebyshev_polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<chebyshev_polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(chebyshev_polynomial<T>(m_nvar,m_degree));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<chebyshev_polynomial<T> > v = evaluate_base(other[j],-1.0,1.0);
        for(int i=0; i<=m_degree; i++)
            base[j][i] = v[i];
    }

    //composing
    chebyshev_polynomial<T> res(m_nvar,m_degree);
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = this->get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            chebyshev_polynomial<T> prod(m_nvar,m_degree);
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

template <class T>
std::vector<T> chebyshev_polynomial<T>::cheb_approximation(T (*f)(T x), const T a, const T b){
    int n = chebyshev_polynomial<T>::MAX_DEGREE;
    std::vector<T> res(n+1), d(n+1);
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

    fac = 2.0/(n+1);

    for (int j = 0; j <= n; j++)
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

// DIRECT MULTIPLICATION
template <class T>
chebyshev_polynomial<T> chebyshev_polynomial<T> :: direct_multiplication(const chebyshev_polynomial<T> &x0, const chebyshev_polynomial<T> &x1){
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

// private routine for direct multiplication
template <class T>
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
