/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
----------------- e-mail:  annalisa.riccardi@strath.ac.uk -------------------
----------------------- Author:  Annalisa Riccardi --------------------------
*/

#include "Sampling/base_sampling.h"

using namespace smart;
using namespace sampling;

template <class T>
base_sampling<T>::base_sampling(const unsigned int &dim, const string &name): m_name(name), m_dim(dim){
}

template <class T>
base_sampling<T>::~base_sampling(){

}

template class base_sampling<double>;
template class base_sampling<float>;
template class base_sampling<long double>;
