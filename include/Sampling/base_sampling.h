#ifndef BASE_SAMPLING_H
#define BASE_SAMPLING_H


/******************************************************************************
 *                            BASE_SAMPLING_H                                 *
 *             Base sampling abstract class of the SMART-UQ toolbox           *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
----------------- e-mail: carlos.ortega@strath.ac.uk ------------------------
----------------------- Author: Carlos Ortega Absil -------------------------
*/

#include <vector>
#include "exception.h"

using namespace std;

namespace smart{
namespace sampling{


    template <class T>
    class base_sampling
    {
        public:

            base_sampling(const unsigned int &dim, const string &name);

            virtual ~base_sampling();

            virtual std::vector<T> operator()() const = 0;

            virtual std::vector<T> operator()(const unsigned int &n) const = 0;

        protected:
            string m_name;
            unsigned int m_dim; 

    };

}}


#endif // BASE_SAMPLING_H
