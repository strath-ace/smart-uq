// # include <cstdlib>
// # include <cmath>
// # include <iostream>
// # include <fstream>
// # include <iomanip>

# include "Non-Intrusive/utils.h"
# include "Non-Intrusive/sparse_grid_index.h"
# include "Non-Intrusive/sparse_grid_dataset.h"
# include "Non-Intrusive/chebyshev_polynomial.h"
# include "Non-Intrusive/multivariate_polynomials.h"

#include <iostream>
#include <Eigen/Dense>

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

using namespace std;
using namespace Eigen;

//****************************************************************************80

int main ()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MULTIVARIATE_CHEBYSHEV_POLYNOMIAL.
//
//  Discussion:
//
//    MAIN tests the MULTIVARIATE_CHEBYSHEV_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{
  timestamp ( );
  cout << "\n";
  cout << "TCHEBYCHEff_POLYNOMIAL\n";
  cout << "  C++ version\n";
  cout << "  Test the CHEBYSHEV_POLYNOMIAL library.\n";

  //  test01 ( );
  
  //  test02( );

  // test03 ( );

  // test04 ( );

  test05( );

  cout << "\n";
  cout << "TCHEBYCHEff_POLYNOMIAL\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 (  )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests the SPASE_GRID_INDEX module.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{ 
  int dim_num;
  int level_max;
  vector< vector<int> > idx;

  cout << "\nInsert spatial dimension:  " ;
  cin >> dim_num;

  cout << "\nInsert level of the approximation:  " ;
  cin >> level_max;

  // Compute multi-index for the sparse grid
  idx = sparse_grid_index( dim_num, level_max);

  // Write on screen
  cout << "\n List of indece:" << endl;
  print_rows(idx, "");

  cout << "\n Number of sequences: " << idx.size() << endl;

}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_CC_DATASET.
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

  int dim_num;
  int level_max;
  vector< vector<int> > idx;

  cout << "\nInsert spatial dimension:  " ;
  cin >> dim_num;

  cout << "\nInsert level of the approximation:  " ;
  cin >> level_max;

  int level_min;
  int point_num;
  string r_filename;
  string w_filename;
  //  double weight_sum;
  string x_filename;
  double *r, *w, *x;

  //  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  //  cout << "  associated with a sparse grid derived from a Smolyak\n";
  //  cout << "  construction based on 1D Clenshaw-Curtis rules.\n";

  //  cout << "    DIM_NUM, the spatial dimension.\n";
  //  cout << "    (typically in the range of 2 to 10)\n";

  //  cout << "    LEVEL_MAX, the level of the sparse grid.\n";
  //  cout << "    (typically in the range of 0, 1, 2, 3, ...\n";

  // cout << "  Output from the program includes:\n";
  // cout << "    * A printed table of the abscissas and weights.\n";
  // cout << "    * A set of 3 files that define the quadrature rule.\n";
  // cout << "    (1) cc_d?_level?_r.txt, the ranges;\n";
  // cout << "    (2) cc_d?_level?_w.txt, the weights;\n";
  // cout << "    (3) cc_d?_level?_x.txt, the abscissas.\n";  

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "  LEVEL_MIN is = " << level_min << "\n";
  cout << "  LEVEL_MAX is = " << level_max << "\n";
// 
//  How many distinct points will there be?
//
  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  The number of distinct abscissas in the\n";
  cout << "  quadrature rule is determined from the spatial\n";
  cout << "  dimension DIM_NUM and the level LEVEL_MAX.\n";
  cout << "  For the given input, this value will be = " << point_num << "\n";
//
//  Allocate memory.
//
  r = new double[dim_num*2];
  w = new double[point_num];
  x = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  for ( int dim = 0; dim < dim_num; dim++ )
  {
    r[dim+0*dim_num] = -1.0;
    r[dim+1*dim_num] = +1.0;
  }

  sparse_grid_cc ( dim_num, level_max, point_num, w, x );

  // r8mat_transpose_print_some ( dim_num, point_num, x, 1, 1, dim_num, 
  //   10, "  First 10 grid points:" );

  // r8vec_print_some ( point_num, w, 1, 10, "  First 10 grid weights:" );

  // weight_sum = 0.0;
  // for ( point = 0; point < point_num; point++ )
  // {
  //   weight_sum = weight_sum + w[point];
  // }

  // cout << "\n";
  // cout << "  Weights sum to   " << weight_sum << "\n";
  // cout << "  Correct value is " << pow ( 2.0, dim_num ) << "\n";
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

  cout << "\n";
  cout << "SPARSE_GRID_CC_DATASET\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

}
//****************************************************************************80

void test03 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests the MULTIVARIATE_TCHEBYCHEFF_MODULE (sparse basis)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    9 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{
  sparse_multivariate_polynomials( dim_num, level_max );
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests the EIGEN library
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
  MatrixXf A = MatrixXf::Random(3, 2);
  cout << "Here is the matrix A:\n" << A << endl;
  VectorXf b = VectorXf::Random(3);
  VectorXf x;
  x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
  cout << "Here is the right hand side b:\n" << b << endl;
  cout << "The least-squares solution is:\n" << x << endl;
  cout << endl << A*x << endl;
  cout << endl << A*x-b << endl;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 test the MULTIVARIATE_TCHEBYCHEFF_MODULE (fill basis)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 October 2015
//
//  Author:
//
//    Chiara Tardioli
//
{
  
  int dim_num;
  int deg_max;
  int point_num;

  cout << "\nInsert spatial dimension:  " ;
  cin >> dim_num;

  cout << "\nInsert degree of expansion: ";
  cin >> deg_max;

  cout << "\nInsert number of points:  " ;
  cin >> point_num;

  vector< vector<double> > pts;
  //  double *x;
  for ( int j=0; j<dim_num; j++)
    {
      pts.push_back( vector<double>(point_num) );
      //      x = r8vec_linspace_new ( point_num, -1., 1. ) ;
      for ( int i=0; i<point_num; i++)
	  pts[j][i] = ((double) rand() / (RAND_MAX)) * 2.0 - 1.0;
    }
  
  vector<double> trues(point_num);
  double sum;
  for ( int i=0; i<point_num; i++)
    {
      sum = 0;
      for ( int j=0; j<dim_num; j++)
	//	  sum += 2.0 * pts[j][i]*pts[j][i] - 1.0; //pow(x[i*dim_num+j],2);
	sum += pts[j][i]*pts[j][i];
      trues[i] = cos(sum);
    }

  vector<double> coe;
  coe = full_multivariate_polynomials ( dim_num, deg_max, pts, trues );

  cout << "----------------------" << endl;
  
  // Print on screen
  //  print_double_matrix(pts, "", "Sample points");
  //  print_double_vector(trues, "", "Trues values");
  //  print_double_vector(coe, "", "Coefficients");

  // Save on file.
  print_double_matrix(pts, "pts.txt", "Sample points");
  print_double_vector(trues, "trues.txt", "Trues values");
  print_double_vector(coe, "coe.txt", "Coefficients");
}
