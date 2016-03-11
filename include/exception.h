/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef SMART_EXCEPTIONS_H
#define SMART_EXCEPTIONS_H

#include <exception>
#include <cassert>
#include <iostream>
#include <string>

#define _EXCEPTION_QUOTEME(x) #x
#define EXCEPTION_QUOTEME(x) _EXCEPTION_QUOTEME(x)
#define EXCEPTION_EXCTOR(s) ((std::string(__FILE__ "," EXCEPTION_QUOTEME(__LINE__) ": ") + s) + ".")
#define EX_THROW(s) (throw smart_exception(EXCEPTION_EXCTOR(s)))

#define smart_throw(s) EX_THROW(s)

#ifdef NDEBUG
        #define smart_assert(condition){ \
            if(!(condition)) { \
                std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__;
                std::cerr << " inside " << __FUNCTION__ << std::endl;
                std::cerr << "Condition: " << #condition;
                abort();
            } \
        }
#else
        #define smart_assert(condition) do { \
            if(!(condition)) { \
                smart_exception("assertion error"); \
            } \
        } while(0)
#endif


class smart_exception: public std::exception {
	public:
		smart_exception(const std::string &s):m_what(s) {}
		virtual const char *what() const throw() {
			return m_what.c_str();
		}
		virtual ~smart_exception() throw() {}
	protected:
		std::string m_what;
};

#endif
