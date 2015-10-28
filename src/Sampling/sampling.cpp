/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
----------------- e-mail:  carlos.ortega@strath.ac.uk -----------------------
----------------------- Author:  Carlos Ortega Absil ------------------------
*/
#include "Sampling/sampling.h"
#include "Sampling/initializers.h"
#include <time.h>

using namespace smart;
using namespace sampling;

/// SOBOL Constructor
/**
 * @param[in] dim dimension of the hypercube
 * @param[in] count starting point of the sequence. choosing 0 wil add the point x=0
 * @throws value_error if dim not in [1,1111]
*/
template <class T>
sobol<T>::sobol(unsigned int dim, unsigned int count) : m_dim(dim), m_count(count), m_dim_num_save(0), m_initialized(false), m_maxcol(62), m_seed_save(-1), recipd(0), lastq(), poly(), v(){
    if (dim >1111 || dim <1) {
      std::cout << "This Sobol sequence can have dimensions [1,1111]" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

/// Operator ()
/**
 * Returns the next point in the sequence
 *
 * @return an std::vector<T> containing the next point
 */

template <class T>
std::vector<T> sobol<T>::operator()() {
  std::vector<T> retval(m_dim,0.0);
  signed long long int seed= (signed long long int) m_count;
  i8_sobol(m_dim, &seed, &retval[0]);
  m_count= (unsigned short int) seed;
  return retval;
}
/// Operator (unsigned int n)
/**
 * Returns the n-th point in the sequence
 *
 * @param[in] n the point along the sequence to be returned
 * @return an std::vector<T> containing the n-th point
 */

template <class T>
std::vector<T> sobol<T>::operator()(unsigned int n) {
  m_count = n;
  signed long long int seed= (signed long long int) m_count;
  std::vector<T> retval(m_dim,0.0);
  i8_sobol(m_dim, &seed, &retval[0]);
  m_count= (unsigned short int) seed;
  return retval;
}


/// LHS Constructor
/**
 * @param[in] dim dimension of the hypercube
 * @param[in] number of points to sample
*/

template <class T>
lhs<T>::lhs(unsigned int dim, unsigned int npoints) :  m_dim(dim), m_npoints(npoints), m_initialised(false), m_set(), m_next(0) {
  srand(time(NULL));
  }
/// Operator ()
/**
 * Returns the next point in the sequence
 *
 * @return an std::vector<T> containing the next point
 */

template <class T>
std::vector<T> lhs<T>::operator()() {
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
/**
 * Returns the n-th point in the sequence
 *
 * @param[in] n the point along the sequence to be returned
 * @return an std::vector<T> containing the n-th point
 */

template <class T>
std::vector<T> lhs<T>::operator()(unsigned int n) {
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

//****************************************************************************
//****************************************************************************
//PRIVATE ROUTINES

//****************************************************************************80
template <class T>
int sobol<T>::i8_bit_lo0 ( long long int n )

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
void sobol<T>::i8_sobol ( unsigned int dim_num, long long int *seed, T quasi[ ] )

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
    long long int polyb[DIM_MAX2];
    initialize_polyb(polyb);

    for ( i = 0; i < DIM_MAX2; i++ )
    {
      poly[i]=polyb[i];
      v[i][0]=1;
      for ( j = 1; j < LOG_MAX; j++ )
      {
        v[i][j] = 0;
      }
    }
//
//  Initialize (part of) V
    initialize_v(v);
//  Check parameters.
//
    if ( dim_num < 1 || DIM_MAX2 < dim_num )
    {
      cout << "\n";
      cout << "I8_SOBOL - Fatal error!\n";
      cout << "  The spatial dimension DIM_NUM should satisfy:\n";
      cout << "    1 <= DIM_NUM <= " << DIM_MAX2 << "\n";
      cout << "  But this input value is DIM_NUM = " << dim_num << "\n";
      exit ( 1 );
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
    cout << "\n";
    cout << "I8_SOBOL - Fatal error!\n";
    cout << "  The value of SEED seems to be too large (too many sample points requested!)\n";
    cout << "  SEED =   " << *seed  << "\n";
    cout << "  MAXCOL = " << m_maxcol << "\n";
    cout << "  L =      " << l << "\n";
    exit ( 2 );
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

//****************************************************************************80

template <class T>
unsigned int *lhs<T>::perm_uniform ( unsigned int n)

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
std::vector<T> lhs<T>::latin_random ( unsigned int dim_num, unsigned int point_num)

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

template class sobol<double>;
template class sobol<float>;
template class sobol<long double>;

template class lhs<double>;
template class lhs<float>;
template class lhs<long double>;