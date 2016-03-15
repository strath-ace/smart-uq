/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Sampling/random_sampling.h"

using namespace smart;
using namespace sampling;

/// random Constructor
template <class T>
random_sampling<T>::random_sampling(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b) :  base_sampling<T>(dim,a,b,"Random Sampling"){
  srand(time(NULL));
  }

/// random Deconstructor
template <class T>
random_sampling<T>::~random_sampling(){
}

/// Operator ()
template <class T>
std::vector<T> random_sampling<T>::operator()() const{
  std::vector<T> retval(m_dim,0.0);
  for (int i=0;i<m_dim;i++){
    retval[i]=((T) rand()) / ((T) RAND_MAX + 1.0);
  }
  return this->map(retval);
}

/// Operator (unsigned int n)
template <class T>
std::vector<T> random_sampling<T>::operator()(const unsigned int &n) const{
  return operator()();
}

template class random_sampling<double>;
template class random_sampling<float>;
template class random_sampling<long double>;
