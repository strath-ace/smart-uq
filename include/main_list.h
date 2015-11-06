#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>

#include "Polynomial/chebyshev.h"
#include "Polynomial/chebyshev_functions.h"

#include "Polynomial/newton.h"

#include "Polynomial/canonical.h"
#include "Polynomial/canonical_functions.h"


#include "integrators_chebyshev.h"
#include "f_chebyshev.h"
#include "integrators_canonical.h"
#include "f_canonical.h"

#include <iterator>
#include <algorithm>

#include <sstream>


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
void main_AKP_i_chebyshev();
void main_AKP_i_canonical();
void main_collision_avoidance();
void main_vanderpol();
void main_test_dct();
void main_test_canonical_multiplication();
void main_test_sampling();
void main_tests();
