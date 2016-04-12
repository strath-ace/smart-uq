/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Sampling/base_sampling.h"

using namespace smartuq;
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
