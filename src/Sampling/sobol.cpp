/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Sampling/sobol.h"

using namespace std;
using namespace smartuq;
using namespace sampling;

/// SOBOL Constructor
template <class T>
sobol<T>::sobol(const unsigned int &dim, const std::vector<T>& a, const std::vector<T>& b, const unsigned int &count) : base_sampling<T>(dim,a,b,"Sobol sampling"), m_count(count), m_dim_num_save(0), m_initialized(false), m_maxcol(62), m_seed_save(-1), recipd(0), lastq(), poly(), v(){
    if (dim >1111 || dim <1)
        smart_throw(m_name+": Sobol sequence can have dimensions [1,1111]");
  }

/// SOBOL Deconstructor
template <class T>
sobol<T>::~sobol(){

}

/// Operator ()

template <class T>
std::vector<T> sobol<T>::operator()() const{
  std::vector<T> retval(m_dim,0.0);
  signed long long int seed= (signed long long int) m_count;
  i8_sobol(m_dim, &seed, &retval[0]);
  m_count= (unsigned short int) seed;
  return this->map(retval);
}

/// Operator (unsigned int n)

template <class T>
std::vector<T> sobol<T>::operator()(const unsigned int &n) const{
  m_count = n;
  signed long long int seed= (signed long long int) m_count;
  std::vector<T> retval(m_dim,0.0);
  i8_sobol(m_dim, &seed, &retval[0]);
  m_count= (unsigned short int) seed;
  return this->map(retval);
}


//****************************************************************************
//****************************************************************************
//PRIVATE ROUTINES

//****************************************************************************80
template <class T>
int sobol<T>::i8_bit_lo0 (long long int n ) const

//****************************************************************************80
//
//  Purpose:
//
//    I8_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1 ??
//    1024 10000000000     1
//    1025 10000000001     1 ??
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long long int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I8_BIT_LO0, the position of the low 1 bit.
//
{
  int bit;
  long long int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
//****************************************************************************80

template <class T>
void sobol<T>::i8_sobol (unsigned int dim_num, long long int *seed, T quasi[ ] ) const

//****************************************************************************80
//
//  Purpose:
//
//    I8_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//    This routine uses LONG LONG INT for integers and T for real values.
//
//    Thanks to Steffan Berridge for supplying (twice) the properly
//    formatted V data needed to extend the original routine's dimension
//    limit from 40 to 1111, 05 June 2007.
//
//    Thanks to Francis Dalaudier for pointing out that the range of allowed
//    values of DIM_NUM should start at 1, not 2!  17 February 2009.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2009
//
//  Author:
//
//    FORTRAN77 original version by Bennett Fox.
//    C++ version by John Burkardt
//
//  Reference:
//
//    IA Antonov, VM Saleev,
//    An Economic Method of Computing LP Tau-Sequences,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley, Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Stephen Joe, Frances Kuo
//    Remark on Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 29, Number 1, pages 49-57, March 2003.
//
//    Ilya Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    Ilya Sobol, YL Levitan,
//    The Production of Points Uniformly Distributed in a Multidimensional
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR,
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 1 <= DIM_NUM <= 1111.
//
//    Input/output, long long int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, T QUASI[DIM_NUM], the next quasirandom vector.
//
{
# define DIM_MAX 40
# define DIM_MAX2 1111
# define LOG_MAX 62
//
//  Here, we have commented out the definition of ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so as long as we set MAXCOL (below) to what we expect it should be, we
//  may be able to get around this difficulty.
//  JVB, 24 January 2006.
//
//static long long int atmost = 4611686018427387903;
//
  long long int i;
  bool includ[LOG_MAX];
  long long int j;
  long long int j2;
  long long int k;
  long long int l;

  long long int m;
  long long int newv;


  long long int seed_temp;

 if ( !m_initialized || dim_num != m_dim_num_save )
 {
    m_initialized = true;

    for (i = 0; i < DIM_MAX2; i++ )
    {
      poly[i]=constants::smart_polyb[i];
      for (j = 1; j < LOG_MAX; j++ )
      {
        v[i][j] = constants::smart_v[i][j];
      }
    }

//  Check parameters.
//
    if ( dim_num < 1 || DIM_MAX2 < dim_num )
    {
      std::stringstream message;
      message << "\n";
      message << "I8_SOBOL - Fatal error!\n";
      message << "  The spatial dimension DIM_NUM should satisfy:\n";
      message << "    1 <= DIM_NUM <= " << DIM_MAX2 << "\n";
      message << "  But this input value is DIM_NUM = " << dim_num << "\n";
      smart_throw(m_name+message.str());
    }

    m_dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
//  Here, we have short-circuited the computation of MAXCOL from ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so if we know what the answer should be we can try to get it this way!
//  JVB, 24 January 2006.
//
//  maxcol = i8_bit_hi1 ( atmost );
//
    m_maxcol = 62;
//
//  Initialize row 1 of V.
//
    for ( j = 0; j < m_maxcol; j++ )
    {
      v[0][j] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; 0 <= k; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < m_maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }
        }
        v[i][j] = newv;
      }
    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = m_maxcol - 2; 0 <= j; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( T ) ( 2 * l ) );

  }//end if(!m_initialised)

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == m_seed_save + 1 )
  {
    l = i8_bit_lo0 ( *seed );
  }
  else if ( *seed <= m_seed_save )
  {
    m_seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = m_seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i8_bit_lo0 ( *seed );
  }
  else if ( m_seed_save+1 < *seed )
  {
    for ( seed_temp = m_seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i8_bit_lo0 ( *seed );
  }
//
//  Check that the user is not calling too many times!
//
  if ( m_maxcol < l )
  {
    std::stringstream message;
    message << "\n";
    message << "I8_SOBOL - Fatal error!\n";
    message << "  The value of SEED seems to be too large (too many sample points requested!)\n";
    message << "  SEED =   " << *seed  << "\n";
    message << "  MAXCOL = " << m_maxcol << "\n";
    message << "  L =      " << l << "\n";
    smart_throw(m_name+message.str());
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( T ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  m_seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef DIM_MAX
# undef DIM_MAX2
# undef LOG_MAX
}



//****************************************************************************
//****************************************************************************

template class sobol<double>;
template class sobol<float>;
template class sobol<long double>;

