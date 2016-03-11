/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "Sampling/lhs.h"

using namespace smart;
using namespace sampling;

/// LHS Constructor
template <class T>
lhs<T>::lhs(const unsigned int &dim, const unsigned int &npoints) :  base_sampling<T>(dim,"Latin Hypercube Sampling"), m_npoints(npoints), m_initialised(false), m_set(), m_next(0) {
  srand(time(NULL));
  }

/// LHS Deconstructor
template <class T>
lhs<T>::~lhs(){

}

/// Operator ()
template <class T>
std::vector<T> lhs<T>::operator()() const{
  std::vector<T> retval(m_dim,0.0);
  if (!m_initialised){
    m_set=latin_random(m_dim,m_npoints);
    m_initialised=true;
    m_next=0;
  }
  for (int i=0;i<m_dim;i++){
    retval[i]=m_set[m_next+i*m_npoints];
  }
  m_next++;
  return retval;
}


/// Operator (unsigned int n)
template <class T>
std::vector<T> lhs<T>::operator()(const unsigned int &n) const{
  std::vector<T> retval(m_dim,0.0);
  m_next = n;
  if (!m_initialised){
    m_set=latin_random(m_dim,m_npoints);
    m_initialised=true;
  }
  for (int i=0;i<m_dim;i++){
    retval[i]=m_set[m_next+i*m_npoints];
  }
  m_next++;
  return retval;
}


//PRIVATE ROUTINES
//****************************************************************************80

template <class T>
unsigned int *lhs<T>::perm_uniform (const unsigned int &n) const

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Output, int PERM_UNIFORM[N], a permutation of (BASE, BASE+1, ..., BASE+N-1).
//
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int *p;
  T r;

  p = new unsigned int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n; i++ )
  {
    r= ((T) rand()) / ((T) RAND_MAX + 1.0);
    j = i+ ( unsigned int ) ( r * (T) (n-i));
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }
  return p;
}

//****************************************************************************80

template <class T>
std::vector<T> lhs<T>::latin_random (const unsigned int &dim_num, const unsigned int &point_num) const

//****************************************************************************80
//
//  Purpose:
//
//    LATIN_RANDOM returns points in a Latin Random square.
//
//  Discussion:
//
//    In each spatial dimension, there will be exactly one
//    point whose coordinate value lies between consecutive
//    values in the list:
//
//      ( 0, 1, 2, ..., point_num ) / point_num
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, T X[DIM_NUM,POINT_NUM], the points.
//
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int *perm;
  T r;
  std::vector<T> x(dim_num*point_num,0.0);
//
//  For spatial dimension I,
//    pick a random permutation of 1 to POINT_NUM,
//    force the corresponding I-th components of X to lie in the
//    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
//
  k = 0;
  for ( i = 0; i < dim_num; i++ )
  {
    perm= perm_uniform ( point_num );

    for ( j = 0; j < point_num; j++ )
    {
      r= ((T) rand()) / ((T) RAND_MAX + 1.0);
      x[k] = ( ( ( T ) perm[j] ) + r ) / ( ( T ) point_num );
      k = k + 1;
    }
    delete [] perm;
  }
  return x;
}

//****************************************************************************
//****************************************************************************


template class lhs<double>;
template class lhs<float>;
template class lhs<long double>;
