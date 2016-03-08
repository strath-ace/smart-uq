/******************************************************************************
 *                       SMART EXCEPTIONS                                     *
 ******************************************************************************/

#ifndef SMART_EXCEPTIONS_H
#define SMART_EXCEPTIONS_H

#include <exception>
#include <cassert>
#include <iostream>
#include <string>

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
