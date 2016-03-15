/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef WRAPPER_H
#define WRAPPER_H
#include "config.h"

#ifdef CHEBYSHEV_DCT_MULTIPLICATION
#include <fftw3.h>
#endif

// RELATED TO THE DCT-BASED MULTIPLICATION
#ifdef CHEBYSHEV_DCT_MULTIPLICATION
// kinda wrap fftw_malloc functions with vanilla template specialization
template <class T, size_t n>
struct dct_malloc_impl;
template <class T>
struct dct_malloc_impl<T,sizeof(float)>{
    void operator()(T*& pointer, int length) const { 
        pointer = (T*) fftwf_malloc( sizeof(T) * length);
        for (int i=0;i<length;i++){
            pointer[i]=0;
        }
    }
};
template <class T>
struct dct_malloc_impl<T,sizeof(double)>{
    void operator()(T*& pointer, int length) const { 
        pointer = (T*) fftw_malloc( sizeof(T) * length);
        for (int i=0;i<length;i++){
            pointer[i]=0;
        }
    }
};
template <class T>
struct dct_malloc_impl<T,sizeof(long double)>{
    void operator()(T*& pointer, int length) const { 
        pointer = (T*) fftwl_malloc( sizeof(T) * length);
        for (int i=0;i<length;i++){
            pointer[i]=0;
        }
    }
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
#endif //CHEBYSHEV_DCT_MULTIPLICATION


#endif // WRAPPERS_H
