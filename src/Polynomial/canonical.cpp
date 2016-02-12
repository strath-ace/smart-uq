/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#include "Polynomial/canonical.h"
#include "Eigen/Dense"

using namespace smart;
using namespace polynomial;

template < class T >
canonical_polynomial<T>::canonical_polynomial(const int &vars, const int &order) : base_polynomial<T>("Canonical Polynomial", vars,order){

}

//initialize a 1 degree univariate chebyshev polynomial of the corresponding variable [x1,x2,...]
template < class T >
canonical_polynomial<T>::canonical_polynomial(const int &vars, const int &order, const int &i)  : base_polynomial<T>("Canonical Polynomial", vars,order,i){

}

//initialize a chebyshev polynomial with only the constant term
template < class T >
canonical_polynomial<T>::canonical_polynomial(const int &vars, const int &order, const T &value) : base_polynomial<T>("Canonical Polynomial", vars,order,value){

}


template < class T >
std::string canonical_polynomial<T>::get_basis_name() const
{
    return "";
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator+(const canonical_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }

    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    canonical_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] + other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator-(const canonical_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    int n = this->get_coeffs().size();

    std::vector<T> other_coeffs = other.get_coeffs();
    std::vector<T> coeffs(n);
    canonical_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<n; i++)
        coeffs[i] = m_coeffs[i] - other_coeffs[i];

    res.set_coeffs(coeffs);
    return res;
}

//OPERATOR* OVERLOADING FOR DIRECT MULTIPLICATION
//initialisation of static terms for multiplication
template <class T>
std::vector<int> canonical_polynomial<T>::m_M(1,0);
template <class T>
int canonical_polynomial<T>::m_Mnvar = 0;
template <class T>
int canonical_polynomial<T>::m_Mdegree = 0;
template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator*(const canonical_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    // use M instead of searching index for faster multiplication
    bool use_M = false;
    if(m_nvar == m_Mnvar && m_degree == m_Mdegree) use_M = true;

    canonical_polynomial<T> res(m_nvar,m_degree);
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
canonical_polynomial<T> canonical_polynomial<T>::inv() const{

    canonical_polynomial<T> res(m_nvar,m_degree);
    std::vector<T> range = canonical_polynomial<T>::get_range();
    T a = range[0];
    T b = range[1];

    std::vector<T> inv = approximation_1d(inverse,a,b,m_degree);
    //univariate composition
    std::vector<canonical_polynomial<T> > base = canonical_polynomial<T>::evaluate_base(a,b);
    for (int i=0; i<=m_degree; i++){
        res += base[i]*inv[i];
    }
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator/(const canonical_polynomial<T> &other) const{

    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    canonical_polynomial<T> res(m_nvar,m_degree);
    res = other.inv();

    return res*(*this);
}


template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator+(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] += other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator-(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);

    coeffs[0] -= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator*(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] *= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator/(const T& other) const{

    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);

    for(int i=0; i<coeffs.size(); i++)
        coeffs[i] /= other;

    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator+() const{

    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::operator-() const{
    
    std::vector<T> coeffs=this->get_coeffs();
    canonical_polynomial<T> res(m_nvar,m_degree);
    for (int i=0;i<coeffs.size();i++){
        coeffs[i] = -coeffs[i];
    }
    res.set_coeffs(coeffs);
    return res;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator=(const canonical_polynomial<T> &other){

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
canonical_polynomial<T>& canonical_polynomial<T>::operator=(const T &other){

    std::vector<T> coeffs(m_coeffs.size());
    coeffs[0] = other;

    m_coeffs = coeffs;
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator+=(const canonical_polynomial<T> &other){
    *this = canonical_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator-=(const canonical_polynomial<T> &other){
    *this = canonical_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator*=(const canonical_polynomial<T> &other){
    *this = canonical_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator/=(const canonical_polynomial<T> &other){
    *this = canonical_polynomial<T>::operator/(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator+=(const T& other){
    *this = canonical_polynomial<T>::operator+(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator-=(const T& other){
    *this = canonical_polynomial<T>::operator-(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator*=(const T& other){
    *this = canonical_polynomial<T>::operator*(other);
    return *this;
}

template <class T>
canonical_polynomial<T>& canonical_polynomial<T>::operator/=(const T& other){
    *this = canonical_polynomial<T>::operator/(other);
    return *this;
}

template <class T>
bool canonical_polynomial<T>::operator==(const canonical_polynomial<T> &other) const{

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
bool canonical_polynomial<T>::operator!=(const canonical_polynomial<T> &other) const{
    if(m_nvar!=other.get_nvar()){
        smart_exception(m_name+"Polynomials don't have the same number of variables. They don't belong to the same Algebra");
    }
    if(m_degree!=other.get_degree()){
        smart_exception(m_name+"Polynomials don't have the same order. They don't belong to the same Algebra");
    }

    if(canonical_polynomial<T>::operator==(other)) return false;
    else return true;

}

//1-d Evaluation method
template <class T>
T canonical_polynomial<T>::evaluate(const T &x) const {
    if(m_nvar>1){
        smart_exception(m_name+"(evaluate) Dimension of point must correspond to number of variables of polynomial.");
    }

    return m_coeffs[0]+horner(x,1);

}

//Multivariate Evaluation method
template <class T>
T canonical_polynomial<T>::evaluate(const std::vector<T> &x) const {
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

//Evaluate canonical base 1, x, x2, x3... in a polynomial. It first map x from [a,b] to [-1,1]
template <class T>
std::vector<canonical_polynomial<T> > canonical_polynomial<T>::evaluate_base(const T &a, const T &b) const{

    std::vector<canonical_polynomial<T> > v;

    for(int i=0; i<=m_degree; i++){
        v.push_back(canonical_polynomial<T>(m_nvar,m_degree));
    }

    //mapping the argument from the domain of composition to [-1,1]
    canonical_polynomial<T> mapped(m_nvar,m_degree);
    if(b==a)
        mapped.set_coeffs(0,a);
    else
        mapped = (2.0*(*this)-(a+b))/(b-a);

    v[0] = 1.0;

    for (int i=1; i<=m_degree; i++){
        v[i] = mapped * v[i-1];
    }
    return v;
}

template <class T>
canonical_polynomial<T> canonical_polynomial<T>::composition(const std::vector<canonical_polynomial<T> > &other) const{
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
    std::vector<std::vector<canonical_polynomial<T> > > base;
    for(int j=0; j<m_nvar; j++){
        std::vector<canonical_polynomial<T> > v;
        for(int i=0; i<=m_degree; i++){
            v.push_back(canonical_polynomial<T>(m_nvar,m_degree));
        }
        base.push_back(v);
    }

    //evaluate all basis
    for(int j=0; j<m_nvar; j++){
        //T range = other[j].get_range();
        std::vector<canonical_polynomial<T> > v = other[j].evaluate_base(-1.0,1.0);
        for(int i=0; i<=m_degree; i++){
            base[j][i] = v[i];
        }
    }

    //composing
    canonical_polynomial<T> res(m_nvar,m_degree);
    int count = 0;
    for(int deg=0; deg<=m_degree; deg++){
        for(int i=0; i<m_J[m_nvar][deg]; i++){
            std::vector<int> row = this->get_row(i,deg); //get for example vector (1 0 0) = x, (0 1 0) = y...
            canonical_polynomial<T> prod(m_nvar,m_degree);
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

//takes a vector of chebyshev coeffs, changes basis and assigns the object to the result
template <class T>
void canonical_polynomial<T>::assign_from_chebyshev(const std::vector<T> cheb_coeffs){ //implementation could be WAY more efficient
    
    canonical_polynomial<T> res(m_nvar,m_degree,(T) 0.0);
    int ncoeffs=res.get_coeffs().size();
    if (cheb_coeffs.size()!=ncoeffs){
        smart_exception(m_name+"Chebyshev coefficients provided must correspond to size of the algebra");
    }

    std::vector <canonical_polynomial <T> > term_vector;
    
    for (int i=0;i<ncoeffs;i++){
        term_vector.push_back(canonical_polynomial<T>(m_nvar,m_degree,(T) cheb_coeffs[i]));
    }

    for (int v=0;v<m_nvar;v++){
        canonical_polynomial<T> base2(m_nvar,m_degree,(T) 1.0);
        canonical_polynomial<T> base1(m_nvar,m_degree,(int) v);
        canonical_polynomial<T> x(m_nvar,m_degree,(int) v);
        canonical_polynomial<T> term(m_nvar,m_degree);
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

    *this=res;
    // return *this;
}

////DIRECT INTERPOLATION WITH MATRIX INVERSION
// template <class T>
// std::vector<T> canonical_polynomial<T>::approximation_1d(T (*f)(T x), const T a, const T b, int degree){
//     // int n = canonical_polynomial<T>::MAX_DEGREE;
//     int n = degree+2;
//     Eigen::MatrixXd base_matrix (n+1,n+1);
//     Eigen::VectorXd y(n+1);

//     std::vector<T> res(n+1);
    
//     T pi = 3.141592653589793;
//     T t, x_mapped;

//     for (int k = 0; k <= n; k++)
//     {
//         t = cos(pi*(k+0.5)/(n+1)); //zeros Ti
//         x_mapped = ((1.0+t)*b + (1.0-t)*a)/2.0; //mapped zeros
//         y(k) = f(x_mapped); //evaluate function
//         base_matrix(k,0)=1.0; //build base_matrix
//         for (int kk=1; kk<=n; kk++){
//             base_matrix(k,kk)=base_matrix(k,kk-1)*t;
//         }
//     }

//     Eigen::VectorXd coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
//     for (int k=0; k<=n; k++){
//         res[k]=coe(k);
//     }

//     return res;
// }


////INTERPOLATION IN CHEBYSHEV BASIS + CHANGE OF BASIS
template <class T>
std::vector<T> canonical_polynomial<T>::approximation_1d(T (*f)(T x), const T a, const T b, const int degree){
    
    int n = canonical_polynomial<T>::MAX_DEGREE;
    // int deg = degree;
    int deg = std::min((int) (degree*1.5+1), n);//RULE OF THUMB
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

    //We translate to canonical basis, we take into acount deg+1 Chebyshev
    //terms but we build a Canonical polynomial of degree+1 terms.
    canonical_polynomial<T> result(1,degree,res[0]);
    canonical_polynomial<T> x(1,degree,(int) 0);
    canonical_polynomial<T> base1(1,degree,(int) 0);
    canonical_polynomial<T> base2(1,degree, (T) 1.0);
    canonical_polynomial<T> base(1,degree);

    result+= res[1]*x;
    for (int i=2;i<=deg;i++){
        base=2.0*x*base1-base2;
        base2=base1;
        base1=base;
        result+=res[i]*base;
    }

    return result.get_coeffs();
}

//private routine for 1d evaluation
template <class T>
T canonical_polynomial<T>::horner(T x, int i) const{
    if (i>m_degree) return 0;
    else return x*(m_coeffs[i]+horner(x,i+1));
}

//initialisation of M
template <class T>
void canonical_polynomial<T>::initialize_M(const int nvar, const int degree){
    
    canonical_polynomial<T> poly(nvar,degree);
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
void canonical_polynomial<T>::delete_M(){
    m_M.resize(1);
    m_Mnvar=0;
    m_Mdegree=0;
}

template class canonical_polynomial<double>;
template class canonical_polynomial<float>;
template class canonical_polynomial<long double>;
