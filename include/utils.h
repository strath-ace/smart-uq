#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>
#include <fftw3.h>

const double ZERO = 1e-15;

template <class T>
T inverse(T x){
    if(fabs(x)<=ZERO){
        std::cout<<"ERROR: Division by zero."<<std::endl;
        throw std::exception();
        //exit(EXIT_FAILURE);
    }
    return 1.0/x;
}

//MATH STUFFS
inline int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline int combination(int n, int k)
{
    int max = std::max(n,k);
    int min = std::min(n,k);
    int res = 1;
    int j = 1;

    while(j<=min){
        res *= max+j;
        j++;
    }

    return res/factorial(min);
}


inline void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, int count){
    if (count < item.size()){
        for (int i = 0; i < values.size(); i++) {
            item[count] = values[i];
            int tmp_count = count + 1;
            rep(res, values, item, tmp_count);
        }
    }else{
        res.push_back(item);
    }
}


inline void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res){
    res.clear();

    std::vector<int> item(k);
    rep(res, values, item, 0);
}

// RELATED TO THE DCT-BASED MULTIPLICATION
// kinda wrap fftw_malloc functions with vanilla template specialization
template <class T, size_t n>
struct dct_malloc_impl;
template <class T>
struct dct_malloc_impl<T,sizeof(float)>{
    void operator()(T*& pointer, int length) const { pointer = (T*) fftwf_malloc( sizeof(T) * length);}
};
template <class T>
struct dct_malloc_impl<T,sizeof(double)>{
    void operator()(T*& pointer, int length) const { pointer = (T*) fftw_malloc( sizeof(T) * length);}
};
template <class T>
struct dct_malloc_impl<T,sizeof(long double)>{
    void operator()(T*& pointer, int length) const { pointer = (T*) fftwl_malloc( sizeof(T) * length);}
};
template <class T>
void dct_malloc(T*& pointer, int length) {dct_malloc_impl<T,sizeof(T)>()(pointer, length);}

// wrap fftw_free functions with vanilla template specializations
template <class T, size_t n>
struct dct_free_impl;
template <class T>
struct dct_free_impl<T,sizeof(float)>{
    void operator()(T*& pointer) const { fftwf_free(pointer);}
};
template <class T>
struct dct_free_impl<T,sizeof(double)>{
    void operator()(T*& pointer) const { fftw_free(pointer);}
};
template <class T>
struct dct_free_impl<T,sizeof(long double)>{
    void operator()(T*& pointer) const { fftwl_free(pointer);}
};
template <class T>
void dct_free(T*& pointer) {dct_free_impl<T,sizeof(T)>()(pointer);}

// function with template specialization that performs everything necessary to obtain a dct and deallocate unnecessary stuff at the end
template <class T, size_t n>
struct dct_do_impl;
template <class T>
struct dct_do_impl<T,sizeof(float)>{
    void operator()(int nvar, int *dct_degree,float*& dct) const {
        fftwf_r2r_kind dct_kind[nvar];
        for (int i=0;i<nvar;i++){
            dct_kind[i]=FFTW_REDFT00;
        }
        fftwf_plan plan = fftwf_plan_r2r(nvar, dct_degree, (float*) dct, (float*) dct, dct_kind, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }
};
template <class T>
struct dct_do_impl<T,sizeof(double)>{
    void operator()(int nvar, int *dct_degree,double*& dct) const {
        fftw_r2r_kind dct_kind[nvar];
        for (int i=0;i<nvar;i++){
            dct_kind[i]=FFTW_REDFT00;
        }
        fftw_plan plan = fftw_plan_r2r(nvar, dct_degree, (double*) dct, (double*) dct, dct_kind, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }
};
template <class T>
struct dct_do_impl<T,sizeof(long double)>{
    void operator()(int nvar, int *dct_degree,long double*& dct) const {
        fftwl_r2r_kind dct_kind[nvar];
        for (int i=0;i<nvar;i++){
            dct_kind[i]=FFTW_REDFT00;
        }
        fftwl_plan plan = fftwl_plan_r2r(nvar, dct_degree, (long double*) dct, (long double*) dct, dct_kind, FFTW_ESTIMATE);
        fftwl_execute(plan);
        fftwl_destroy_plan(plan);
    }
};
template <class T>
void dct_do(int nvar, int* dct_degree, T*& dct) {
    dct_do_impl<T,sizeof(T)>()(nvar,dct_degree,dct);
}


//LAPACK METHOD
//    std::vector<T> coeffs = other.get_coeffs();
//    int n = coeffs.size();
//    int nrhs = 1;
//    double B[n][n];
//    double b[1][n];
//    int lda = n;
//    int ldb = n;
//    int ipiv[n];
//    int info;

//    for(int i=0; i<n; i++){
//        b[0][i] = 0.0;
//        for(int j=0; j<=i; j++){
//            if(j==0){//first row and first diagonal
//                B[i][0] = coeffs[i];
//                B[j][i] =  B[i][j];
//            }
//            else if (i==j){//diagonal
//                if(2*i < n)
//                    B[i][i] = 2.0*coeffs[0]+coeffs[2*i];
//                else
//                    B[i][i] = 2.0*coeffs[0];
//            }
//            else{//lower diagonal
//                if((i+j)<n)
//                    B[i][j] = coeffs[fabs(i-j)]+coeffs[i+j];
//                else
//                    B[i][j] = coeffs[fabs(i-j)];
//                B[j][i] =  B[i][j];
//            }
//        }
//    }

//    b[0][0] = 2.0;

//    dgesv_(&n, &nrhs, &B[0][0], &lda, ipiv, &b[0][0], &ldb, &info);

//    // Check for success
//    if(info == 0)
//    {
//       std::vector<T> res_coeffs(n);
//       for(int i=0; i<n; i++)
//           res_coeffs[i] = b[0][i];
//       res_coeffs[0] /= 2.0;
//       res.set_coeffs(res_coeffs);
//    }
//    else
//    {
//       // Write an error message
//       std::cout << "LAPACK dgesv returned error " << info << "\n";
//       exit(EXIT_FAILURE);
//    }


#endif // UTILS_H
