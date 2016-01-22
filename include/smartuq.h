/******************************************************************************
 *                               SMARTUQ_H                                    *
 *            Recursive inclusion of the SMART-UQ headers                     *
 ******************************************************************************/

#ifndef SMART_SMARTUQ_H
#define SMART_SMARTUQ_H

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

#endif /* SMART_SMARTUQ_H */
