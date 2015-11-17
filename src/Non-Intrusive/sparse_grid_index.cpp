/*
 *  ------------------------------------------------------------------------
 *
 *                              M O D U L E
 *
 *         S P A R S E    G R I D    I N D E X    M A T R I X
 *
 *  ------------------------------------------------------------------------
 *
 * Author:
 *   Chiara Tardioli
 *
 * Revided:
 *   October 9, 2015
 *
 */

# include "Non-Intrusive/sparse_grid_index.h"

//****************************************************************************80

vector<int> sumVector( vector<int> a, vector<int> b)

//****************************************************************************80
//
//  Purpose:
//
//    SUMVECTOR return the sum of two integer vectors
//
//  Discussion
//
//    It works only for integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector<int> A, first vector of integer
//
//    Input, vector<int> B, second vector of integer
//
//    Output, vector<int> SUMVECTOR, the vector sum
//
{
  unsigned int n=a.size();
  vector<int> c(n);

  if( n != b.size() )
    {
      cout << "ERROR : vectors with different size!" << endl;
      cout << "SIZE = "<< n << ' ' << b.size() << endl;
      exit(1);
    }
  
  for( unsigned int i=0; i<n; i++)
    c[i] = a[i] + b[i];

  return c;
}
//****************************************************************************80

vector<int> array2vector(int arr[], int n)

//****************************************************************************80
//
//  Purpose:
//
//    ARRAY2VECTOR transfor a class array to a class vector, and copy the integers.
//
//  Discussion
//
//    It works only for integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int ARR, the list of integers in array class
//
//    Input, int N, the lenght of ARR
//
//    Output, vector<int> ARRAY2VECTOR, the corresoponding list of integer in vector class
//
{
  vector<int> v(n);
  for( int i=0; i<n; i++)
    v[i] = arr[i];
  return v;
}
//****************************************************************************80

void printArray(int arr[], int size) 

//****************************************************************************80
//
//  Purpose:
//
//    PRINTARRAY print an array on screen
//
//  Discussion
//
//    It works only for integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int ARR, the list of integers in array class
//
//    Input, int N, the lenght of ARR
//
{
  for ( int i = 0; i < size; i++ )
    cout << arr[i] << " ";
  cout << endl;
}
//****************************************************************************80

void printVector(vector<int>& v, int level, ofstream& myfile) 

//****************************************************************************80
//
//  Purpose:
//
//    PRINTVECTOR print a vector on file if myfile is open, on screen otherwise.
//
//  Discussion
//
//    It works only for integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector<int>& V, the list of integers in vector class
//
//    Input, int LEVEL, the components of vector to print, from left to right, after that it put zeros
//
//    Input, ofstream& MYFILE, the pointer to the open file
//
{
  vector<int>::iterator it;

  if (myfile.is_open())
    {
      myfile << level << " ";
      for(int i = 0; i < level+1; ++i)
	myfile << v[i] << " ";
      for(unsigned int i = level+1; i < v.size(); ++i)
	myfile << 0 << " ";
      myfile << endl;
    }
  else
    { 
      //cout << level << " ";
      for(int i = 0; i < level; ++i)
	cout << v[i] << " ";
      //      for(unsigned int i = level; i < v.size(); ++i)
      //	cout << 0 << " ";
      cout << endl;

      //      for(it = v.begin(); it != v.end(); ++it)
      //      	cout << *it << " ";
      //      cout << endl;
    }
}
//****************************************************************************80

void print_rows(vector< vector<int> > vv, string filnam)

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_ROWS print a matrix row by row if filname is a non-empty string
//
//  Discussion
//
//    It works only for integers. If filname is a non-empty string, it opens the file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector< vector<int> > VV, the integer matrix
//
//    Input, string FILNAM, name of the output file to open (empty for the standard output)
//
{
  vector<int>::iterator it;

  if( filnam.length() )
    {
      ofstream out;
      const char *fil = filnam.c_str();
      out.open(fil);

      for( unsigned int i=0; i<vv.size(); ++i)
	{
	  out << "    " ;
	  for(it = vv[i].begin(); it != vv[i].end(); ++it)
	    out << *it << " ";
	  out << endl;
	}
      out.close();
    }
  else
    {
      for( unsigned int i=0; i<vv.size(); ++i)
	{
	  for(it = vv[i].begin(); it != vv[i].end(); ++it)
	    cout << *it << " ";
	  cout << endl;
	}
    }
}
//****************************************************************************80

vector< vector<int> > part_recursive(int n, vector<int>& v, int level, 
				     vector< vector<int> > vv, int dim)

//****************************************************************************80
//
//  Purpose:
//
//    PART_RECURSIVE compute the partition of N with DIM dimension by recursion
//
//  Discussion
//
//    Of course, it works only for integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int N, number to be partitioned
//
//    Input, vector<int>& v, auxiliar vector of partitions fot the recursion
//
//    Input, int LEVEL, level of the recursion
//
//    Input, vector< vector<int> > VV, all partitions computed so far in a matrix format
//
//    Input, int DIM, number of bins
//
//    Output, vector< vector<int> > PART_RECURSIVE, all partitions in a matrix format
//
{
    int first; /* first is before last */

    if(n<1) return vv; // never goes here, boh!
    v[level]=n;

    if( level+1 == dim )
      {
	vector<int> vaux(dim);
	for(int i=0; i<dim; ++i)
	  vaux[i] = v[i];
	vv.push_back( vaux );
	return vv;
      }

    first=(level==0) ? 1 : v[level-1];
    
    for(int i=first;i<=n/2;i++)
      {
        v[level]=i; /* replace last */
        vv = part_recursive(n-i, v, level+1, vv, dim );
      }

    return vv;
}
//****************************************************************************80

vector< vector<int> > partitions(int N, int dim)

//****************************************************************************80
//
//  Purpose:
//
//    PARTITIONS return the partition of N with DIM bins.
//
//  Discussion
//
//    This is a wrap function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int N, number to be partitioned
//
//    Input, int DIM, number of the bins
//
//    Output, vector< vector<int> > PARTITIONS, all partitions returned in a matrix format
//
{ 
  // Compute partitions.
  vector<int> v(N);
  vector< vector<int> > vv;
  vv = part_recursive( N, v, 0, vv, dim );

  // Return partitions.
  return vv;
}
//****************************************************************************80

vector< vector<int> > permutations( vector<int> v, int dim, vector< vector<int> > vv )

//****************************************************************************80
//
//  Purpose:
//
//    PERMUTATIONS computes all the permutation without repetitions of a vector
//       of integers by recursion.
//
//  Discussion
//
//    It works only for integers. If DIM > size(V), zeros are added.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector<int> V, vector of integer to be permuted
//
//    Input, int DIM, space dimension
//
//    Input, vector< vector<int> > VV, permutation computed so far (needed for recursion)
//
//    Output, vector< vector<int> > PERMUTATION, all partitions returned in a matrix format
//
{  
  int myints[dim];
  vector<int> w(dim);

  for (int i=0; i<dim; i++)
      myints[i] = v[i];
  sort(myints, myints+dim);

  do {
    //printArray(myints, dim);
    w = array2vector(myints, dim);
    vv.push_back( w );
  } while ( next_permutation(myints, myints+dim) );

  return vv;
}
//****************************************************************************80

int emme(int i)

//****************************************************************************80
//
//  Purpose:
//
//    EMME compute the Smolyak multi-index for the size of the 1D sparse grid
//
//  Discussion
//
//    Algorithm:  m(i) = 2^(i-1) + 1 for i>=2, m(1)=1, m(0)=0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int i, level of the spoarse grid
//
//    Output, int EMME, number of points in the 1D sparse grid
//
{
  if( i == 0 )  return 0;
  if( i == 1 )  return 1;
  return pow(2,i-1) + 1;
}
//****************************************************************************80

vector< vector<int> > smolyak_recursive(int n, int dim, int lower[], int upper[],
					vector<int> v, vector< vector<int> > vv)

//****************************************************************************80
//
//  Purpose:
//
//    SMOLYAK_RECURSIVE returns the multi-index matrix according to the nested Smolyak's rule
//
//  Discussion
//
//    See Judd et al. for the algorithmn.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int N, integer counter for recursion
//
//    Input, int DIM, space dimension
//
//    Inout, int LOWER[], minimum indeces in the doloop
//
//    Input int UPPER[], maximum indeces in the doloop
//
//    Input, vector<int> V, vector of indeces filled in the recursion
//
//    Input  vector< vector<int> > VV, stored matrix of indeces
//
//    Output, vector< vector<int> SMOLYAK_RECURSIVE, new matrix of indeces (old + new)
//
{
  if( n >= dim )
    {
      vv.push_back( v );
      return vv;
    }
  for( int i=lower[n]; i<upper[n]+1; i++)
    {
      v[n] = i;
      vv = smolyak_recursive(n+1, dim, lower, upper, v, vv);
    }  
  return vv;
}
//****************************************************************************80

vector< vector<int> > smolyak_indices(vector<int> v, int dim, vector< vector<int> > vv)

//****************************************************************************80
//
//  Purpose:
//
//    SMOLYAK_INDICES returns the multi-index matrix according to the nested Smolyak's rule
//
//  Discussion
//
//    This is a wrap function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, vector<int> V, list of integer (it represent a permutation vector)
//
//    Input, int DIM, space dimension
//
//    Input  vector< vector<int> > VV, previuos stored matrix of indeces
//
//    Output, vector< vector<int> SMOLYAK_RECURSIVE, new matrix of indeces (old + new)
//
{
  // Set bouns for multi-index array.
  int lower[dim];
  int upper[dim];

  for( int i=0; i<dim; i++)
    {
      lower[i] = emme( v[i]-1 ) + 1;
      upper[i] = emme( v[i] );
    }

  // Compute multi-indices and write them on temporary file.
  vector<int> w(dim);
  vv = smolyak_recursive( 0, dim, lower, upper, w, vv );

  // Return list of indices.
  return vv;
}
//****************************************************************************80

vector< vector<int> > sparse_grid_index( int dim, int mu)
 
//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_INDEX compute all the sparse grid indeces for the given level 
//        of expansion and number of variables.
//
//  Discussion
//
//    See Judd et al. for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
//  Parameters:
//
//    Input, int DIM, space dimension
//
//    Input, int MU, level of the sparse grid
//
//    Output, vector< vector<int> > SPARSE_GRID_INDEX, list of multi-index
//            for the required sparse grid 
//
{
  vector< vector<int> > pt, pm, idx;
  vector<int> v(dim);
  
  for(int N=dim; N<=dim+mu; N++)
    {
      // clean up.
      pm.clear();

      //      cout << "----- PARTITIONS " << endl;
      pt = partitions(N, dim);

      //      cout << " ---- PERMUTATIONS " << endl;
      for( unsigned int i=0; i<pt.size(); i++)
	  pm = permutations( pt[i], dim, pm);

      //      cout << "\n ---- SMOLYAK POLINOMIALS: efficient costruction " << endl;
      for( unsigned int i=0; i<pm.size(); i++)
	  idx = smolyak_indices(pm[i], dim, idx);
    }

  // Put each index in the range [0,2^mu]
  for( unsigned int i=0; i<idx.size(); i++)
      for( int j=0; j<dim; j++)
	       idx[i][j] -= 1;

  // Reorder index (Reference: Giorgilli e Sansottera)
  vector< vector<int> > Jtab, Ntab;
  int deg_max = pow(2,mu);
  Jtab = initialize_J( dim, deg_max);
  Ntab = initialize_N( dim, deg_max, Jtab);

  int index, deg, Jdeg;
  int tot = Ntab[dim][deg_max];
  vector<int> idx_col(tot);
  for ( int i=0; i<tot; i++ )
    idx_col[i] = 0;

  for ( unsigned int i=0; i<idx.size(); i++ )
    {
      deg=0;
      for ( int j=0; j<dim; j++)
	deg += idx[i][j];

      if (deg==0)
	Jdeg = 0;
      else
	Jdeg = Ntab[dim][deg-1];

      index = vector_to_idx( idx[i], dim, Ntab) + Jdeg;
      idx_col[index] = 1;
    }

  int i_prev=0, count=0;
  for(int deg=0; deg<=deg_max; ++deg) {
      for(int j=0; j<Jtab[dim][deg]; ++j)
	{
	  if (deg==0)
	    i_prev = 0;
	  else
	    i_prev = Ntab[dim][deg-1];

	  if ( idx_col[i_prev + j] == 1)
	    {
	      idx[count] = idx_to_vector( j, deg, dim, Ntab);
	      count += 1;
	    }
    }
	}

  // // Write on file.
  // string filename;
  // ostringstream convert1, convert2;   // stream used for the conversion

  // string dimension;
  // convert1 << dim; 
  // dimension = convert1.str();

  // string level;
  // convert2 << mu; 
  // level = convert2.str();

  // filename = "index_" + dimension + "_" + level + ".txt";

  // print_rows( idx, filename );

  return idx;
}
//****************************************************************************80

vector< vector<int> > initialize_J( int dim, int deg)

//****************************************************************************80
//
//  Purpose:
//
//    initialize_J from ALCOR to handle index order.
//
//  Discussion
//
//    See Giorgilli and Sansottera for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   20 September 2014
//
//  Author:
//
//    Marco Sansottera
//
//  Parameters:
//
//    Input, int DIM, space dimension
//
//    Input, int DEG, maximum degree of approximation
//
//    Output, vector< vector<int> > Jtab
//
{
  int i,j;
  vector< vector<int> > Jtab(dim+1, vector<int>(deg+1));
  
  //fill J
  for(j = 0; j <= deg; ++j)
        Jtab[1][j] = 1;
  for(i = 2; i <= dim; ++i) 
    {
      Jtab[i][0] = 1;
      for(j = 1; j <= deg; ++j)
	Jtab[i][j] = Jtab[i][j-1]+Jtab[i-1][j];
    }
    
  return Jtab;
}
//****************************************************************************80

vector< vector<int> > initialize_N( int dim, int deg, vector< vector<int> > Jtab)

//****************************************************************************80
//
//  Purpose:
//
//    initialize_N from ALCOR to handle index order.
//
//  Discussion
//
//    See Giorgilli and Sansottera for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   20 September 2014
//
//  Author:
//
//    Marco Sansottera
//
//  Parameters:
//
//    Input, int DIM, space dimension
//
//    Input, int DEG, maximum degree of approximation
//
//    Output, vector< vector<int> > Ntab
//
{
  int i,j;
  vector< vector<int> > Ntab(dim+1, vector<int>(deg+1));
  
  //fill N
  for(i = 1; i <= dim; ++i)
    Ntab[i][0] = Jtab[i][0];
  for(i = 1; i <= dim; ++i) 
    {
      for(j = 1; j <= deg; ++j)
	Ntab[i][j] = Ntab[i][j-1] + Jtab[i][j];
    }
      
  return Ntab;
}
//****************************************************************************80

int vector_to_idx(vector<int> k, int dim, vector< vector<int> > Ntab)

//****************************************************************************80
//
//  Purpose:
//
//    VECTOR_TO_IDX from ALCOR to handle index order.
//
//  Discussion
//
//    See Giorgilli and Sansottera for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   20 September 2014
//
//  Author:
//
//    Marco Sansottera
//
//  Parameters:
//
//    Input, vector<int> k, vector of integers
//
//    Input, int DIM, space dimension
//
//    Input, vector< vector<int> > Ntab, N matrix
//
//    Output, int IDX, index integer of the homogeneous polynomial
//
{
    int i, j, l;
    int idx = 0;

    for(i = dim-1; (i>=0) && (k[i]==0); --i);

    if(i < 0)
        return idx;
    j = dim-i;
    l = -1;
    while(i >= 0) {
        l += k[i--];
        idx += Ntab[j++][l];
    }
    idx -= Ntab[dim][l];

    return idx;
}
//****************************************************************************80

vector<int> idx_to_vector( int idx, int deg, int dim, vector< vector<int> > Ntab)

//****************************************************************************80
//
//  Purpose:
//
//    IDX_TO_VECTOR from ALCOR to handle index order.
//
//  Discussion
//
//    See Giorgilli and Sansottera for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   20 September 2014
//
//  Author:
//
//    Marco Sansottera
//
//  Parameters:
//
//    Input, int IDX, index integer of the homogeneous polynomial
//
//    Input, int DEG, degree of the subgroup
//
//    Input, int DIM, space dimension
//
//    Input, vector< vector<int> > Ntab, N matrix
//
//    Output, vector<int> k, vector of integers
//
{
    int i, j, l;
    int idx_tmp = idx;
    vector<int> k(dim);

    if(deg > 0)
        idx_tmp += Ntab[dim][deg-1];

    l = dim;
    for(i = 0; i < dim; ++i) {
        if(idx_tmp == 0) {
            for(j = i; j < dim; k[j++] = 0);
            return k;
        }
        if(i == 0)
            j = deg;
        else {
            for(j = 0; idx_tmp >= Ntab[l][j]; ++j);
            k[i-1] -= j;
        }
        k[i] = j;
        idx_tmp -= Ntab[l--][j-1];
    }

    return k;
}
