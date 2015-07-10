#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <math.h>

#include "chebyshev_polynomial.h"
#include "elementary_functions.h"
#include "integrators.h"
#include "f.h"


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


void main_accuracy();
void main_AKP();
void main_collision_avoidance();
void main_vanderpol();

