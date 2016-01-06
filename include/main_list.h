#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <sstream>
//chebyshev algebra
#include "Polynomial/chebyshev.h"
#include "Polynomial/chebyshev_functions.h"
//newton algebra
#include "Polynomial/newton.h"
//canonical algebra
#include "Polynomial/canonical.h"
#include "Polynomial/canonical_functions.h"
//chebyshev testcases
#include "integrators_chebyshev.h"
#include "f_chebyshev.h"
//canonical testcases
#include "integrators_canonical.h"
#include "f_canonical.h"
//non-intrusive testcases
#include "utils.h"
#include "Sampling/sampling.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "f_ni.h"
#include "integrators_ni.h"
//non-intrusive sparse additional routines
# include "Non-Intrusive/utils.h"
# include "Non-Intrusive/sparse_grid_index.h"
# include "Non-Intrusive/sparse_grid_dataset.h"
# include "Non-Intrusive/chebyshev_polynomial.h"
# include "Non-Intrusive/multivariate_polynomials.h"



namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}


void main_multiphase();
void main_test_accuracy();
void main_AKP_ni();
void main_AKP_ni_sparse();
void main_AKP_i_chebyshev();
void main_AKP_i_canonical();
void main_collision_avoidance();
void main_vanderpol_ni();
void main_vanderpol_ni_sparse();
void main_vanderpol_i_chebyshev();
void main_vanderpol_i_canonical();

void main_lotka_volterra_ni();
void main_lotka_volterra_ni_sparse();
void main_lotka_volterra_i_chebyshev();
void main_lotka_volterra_i_canonical();

void main_spring_mass_ni();

void main_test_dct();
void main_test_canonical_multiplication();
void main_test_sampling();
void main_tests();
