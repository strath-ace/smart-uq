#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <math.h>
#include <stdio.h>
# include <iomanip>     // setw
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

# include "Non-Intrusive/multivariate_polynomials.h"
# include "Non-Intrusive/chebyshev_polynomial.h"
# include "Non-Intrusive/sparse_grid_index.h"
# include "Non-Intrusive/sparse_grid_dataset.h"

//****************************************************************************80

vector<double> full_multivariate_polynomials ( int dim_num, int deg_max, 
			 vector< vector<double> > pts, vector<double> trues )

//****************************************************************************80
//
//  Purpose:
//
//    FULL_MULTIVARIATE_POLYNOMIALS solve the Tchebycheff coefficients via 
//          least-square methods.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{
  //
  // Number of sample points.
  // 
  if (pts.size() < 0)
    {
      cout<<"\n Maximum allowed polynomial degree is 20\n";
      exit(EXIT_FAILURE);
    }
  int point_num = pts[0].size();
  cout << "\n No. of sample points = " << point_num << endl;

  //
  // Compute T_(n,x_i^k) for 0<=n<deg_max and 1<=i<dim_num, 1<=k<m.
  //
  double *v, *x;
  vector<double> T(dim_num*point_num*(deg_max+1));
  x = new double[point_num];
  for( int i=0; i<dim_num; i++)
    {
      for ( int k=0; k<point_num; k++ )
	  x[k] = pts[i][k];
      v = t_polynomial ( point_num, deg_max, x );

      for( int j=0; j<point_num*(deg_max+1); j++)
	  T[i*point_num*(deg_max+1)+j] = v[j]; // Add the row to the main vector
    }
  delete [] x;

  //
  // Total number of elements : dim_num*point_num*(deg_max+1).
  //
  cout << " Total number of elements : " << T.size() << endl;

  //
  // Initialize INDEX MATRIX (full)
  //
  vector< vector<int> > Jtab, Ntab;
  Jtab = initialize_J( dim_num, deg_max);
  Ntab = initialize_N( dim_num, deg_max, Jtab);

  vector< vector<int> > idx;
  int count=0;
  for(int i=0; i<=deg_max; ++i) 
      for(int j=0; j<Jtab[dim_num][i]; ++j)
	{
	  idx.push_back( idx_to_vector( j, i, dim_num, Ntab) );
	  count += 1;
	}
  int len=idx.size();
 
  cout << " Degree = " << deg_max << "\n Dimension = " << dim_num << endl;
  cout << " No. of coeffs = " << len << endl;

  //
  // Define matrix H/HP.
  //
  vector< vector<double> >H;
  vector<double> row(len);
  double prod;

  for ( int j=0; j<point_num; ++j )
    {
      for ( int k=0; k<len; ++k )
	{
	  prod=1;
	  for ( int i=0; i<dim_num; ++i)
	      prod *= T[i*point_num*(deg_max+1) + idx[k][i]*point_num + j];
	  row[k] = prod;
	}
      H.push_back( row );
    }

  MatrixXf HP(point_num, len);
  for( int i=0; i<point_num; ++i)
      for( int j=0; j<len; ++j )
	HP(i,j) = H[i][j];

  VectorXf b(point_num);
  for( int i=0; i<point_num; ++i)
    b(i) = trues[i];

  //
  // Invert H matrix and solve for the unknow coefficients
  //
  VectorXf coe,err;
  coe = HP.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
  err = HP*coe-b;
  double rmse = sqrt(err.adjoint()*err);

  //  cout << "Here is the matrix HP:\n" << HP << endl;
  //  cout << "Here is the right hand side b:\n" << b << endl;
  //  cout << "The least-squares solution is:\n" << coe << endl;
  cout << endl << "Root mean square error is:\n" << rmse << endl;

  vector<double> vcoe;
  for( int i=0; i<len; ++i)
    vcoe.push_back( coe[i] );

  return vcoe;
}
//****************************************************************************80

void sparse_multivariate_polynomials ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_MULTIVARIATE_POLYNOMIALS solve the Tchebycheff coefficients via 
//          least-square methods.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{
  // 
  //  How many distinct points will there be?
  //
  int point_num;
  point_num = sparse_grid_cfn_size ( dim_num, level_max );
  cout << "no. pts = " << point_num << endl;

  //
  // Compute sparse grid points.
  //
  double *xpr, *x;
  xpr = sparse_grid_points ( dim_num, level_max, point_num );

  //
  // Compute T_(n,x_i^k) for 0<=n<deg_max and 1<=i<dim_num, 1<=k<m.
  //
  const double deg_max = pow(2,level_max);
  double *v;
  vector<double> T(dim_num*point_num*(deg_max+1));
  x = new double[point_num];
  for( int i=0; i<dim_num; i++)
    {
      for ( int k=0; k<point_num; k++ )
	      x[k] = xpr[k*dim_num+i];
      
      v = t_polynomial ( point_num, deg_max, x );

      for( int j=0; j<point_num*(deg_max+1); j++)
        T[i*point_num*(deg_max+1)+j] = v[j]; // Add the row to the main vector
    }
  delete [] x;

  //
  // Total number of elements : dim_num*point_num*(deg_max+1).
  //
  cout << "Total number of elements : " << T.size() << endl;

  //  
  // Check: print on screen
  //
  // for( int i=0; i<dim_num; ++i)
  //   {
  //     cout << " variable " << i << endl;
  //     for( int j=0; j<point_num; ++j )
  // 	{
  // 	  for( int k=0; k<deg_max+1; ++k )
  // 	    cout << T[i*point_num*(deg_max+1)+k*point_num+j] << " ";
  // 	  cout << endl;
  // 	}
  //   }

  //
  // Define matrix H/HP.
  //
  vector< vector<int> > idx;
  idx = sparse_grid_index( dim_num, level_max);

  int len=idx.size();
  vector< vector<double> >H;
  vector<double> row(len);
  double prod;

  for ( int j=0; j<point_num; ++j ){
      for ( int k=0; k<len; ++k ){
	        prod=1;
	        for ( int i=0; i<dim_num; ++i)
	            prod *= T[i*point_num*(deg_max+1) + idx[k][i]*point_num + j];
	            row[k] = prod;
	    }
      H.push_back( row );
  }

  MatrixXf HP(point_num, len);
  for( int i=0; i<point_num; ++i)
      for( int j=0; j<len; ++j )
	       HP(i,j) = H[i][j];
  

  vector<double> trues(point_num);
  for ( int i=0; i<point_num; i++)
    {
      prod = 1;
      for ( int j=0; j<dim_num; j++)
	prod += pow(xpr[i*dim_num+j],2);
      trues[i] = prod;
    }

  VectorXf b(len);  //HERE: len=point_num
  for( int i=0; i<point_num; ++i)
    b(i) = trues[i];

  //
  // Invert H matrix and solve for the unknow coefficients
  //
  VectorXf coe,err;
  coe = HP.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
  err = HP*coe-b;
  double rmse = sqrt(err.adjoint()*err);

  cout << "Here is the matrix HP:\n" << HP << endl;
  cout << "Here is the right hand side b:\n" << b << endl;
  cout << "The least-squares solution is:\n" << coe << endl;
  cout << endl << "Root mean square error is:\n" << rmse << endl;

  string filnam;
  ofstream out;
  filnam = "coe.txt";
  const char *fil1 = filnam.c_str();
  out.open(fil1);
  out << coe << endl;
  out.close();

  filnam = "trues.txt";
  const char *fil2 = filnam.c_str();
  out.open(fil2);
  out << b   << endl;
  out.close();
}
//****************************************************************************80

double* sparse_grid_points ( int dim_num, int level_max, int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_POINTS is the wrap function for SPARSE_GRID_CC_DATASET.
//
//  Discussion:
//
//    This program computes a sparse grid quadrature rule based on 1D
//    Clenshaw-Curtis rules and writes it to a file.
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
{
  int dim;
  int level_min;
  string r_filename;
  string w_filename;
  //  double weight_sum;
  string x_filename;
  double *r, *w, *x;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "  LEVEL_MIN is = " << level_min << "\n";
  cout << "  LEVEL_MAX is = " << level_max << "\n";
//
//  Allocate memory.
//
  r = new double[dim_num*2];
  w = new double[point_num];
  x = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  for ( dim = 0; dim < dim_num; dim++ )
  {
    r[dim+0*dim_num] = -1.0;
    r[dim+1*dim_num] = +1.0;
  }

  sparse_grid_cc ( dim_num, level_max, point_num, w, x );
//
//  Construct appropriate file names.
//
  r_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
  w_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
  x_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
//
//  Write the rule to files.
//
  cout << "\n";
  cout << "  Creating R file = \"" << r_filename << "\".\n";

  r8mat_write ( r_filename, dim_num, 2, r );

  cout << "  Creating W file = \"" << w_filename << "\".\n";

  r8mat_write ( w_filename, 1, point_num, w );

  cout << "  Creating X file = \"" << x_filename << "\".\n";

  r8mat_write ( x_filename, dim_num, point_num, x );
//
//  Terminate.
//
//  delete [] r;
//  delete [] w;
//  delete [] x;

  return x;
}
