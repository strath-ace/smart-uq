/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include "base_polynomial.h"
#include "config.h"

#ifdef CHEBYSHEV_DCT_MULTIPLICATION
#include <fftw3.h>
#endif

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
    /**
     *@brief Chebyshev polynomial
     *
     * The chebyshev_polynomial class implements the algebra between chebyshev polynomials and additional methods for approximation
     * of given function, composition and polynomial evaluation
     */
    class chebyshev_polynomial: public base_polynomial <T> {
	
	// // for class template inheritance
	private:
	using base_polynomial<T>::m_name;
	using base_polynomial<T>::m_coeffs;
	using base_polynomial<T>::m_degree;
	using base_polynomial<T>::m_nvar;
	using base_polynomial<T>::m_J;
	using base_polynomial<T>::m_N;
    using base_polynomial<T>::m_a;
    using base_polynomial<T>::m_b;
    using base_polynomial<T>::m_monomial_base;

    public:
        /**
         * @brief MAX_DEGREE maximum polynomial degree expansion used for univariate polynomial approximation
         */
        static const int MAX_DEGREE=100;
        /**
         * @brief chebyshev_polynomial constructor
         *
         * The constructor initializes a null polynomial for a given number of variables and degree.
         * It is possible through the constructor to define the range of each variable. If range is not specified [-1,1]
         * is assumed
         * @param vars number of variables
         * @param order polynomial maximum degree
         * @param a vector containing the lower bound for each variable (default -1)
         * @param b vector containing the upper bound for each variable (default 1)
         * @param monomial flag to detect the change of base (default false)
         */
        chebyshev_polynomial(const int &vars, const int &order, const std::vector<T> &a=std::vector<T>(), const std::vector<T> &b=std::vector<T>(), const bool &monomial=false);
        /**
         * @brief chebyshev_polynomial constructor
         *
         * The constructor initializes the first order polynomial degree in the i-th variable.
         * The coefficients for degree up to the specified order are allocated and set to zero.
         * The variable is mapped from [a,b] to [-1,1]
         * @param vars number of variables
         * @param order polynomial maximum degree
         * @param i index of first degree variable
         * @param a variable lower bound (default value -1)
         * @param b variable upper bound (default vale 1)
         * @param monomial flag to detect the change of base (default false)
         */
        chebyshev_polynomial(const int &vars, const int &order, const int &i, const T &a=-1.0, const T &b=1.0, const bool &monomial=false);
        /**
         * @brief chebyshev_polynomial
         *
         * The constructor initialize the zero order polynomial degree with given value.
         * The coefficients for degree up to the specified order are allocated and set to zero
         * @param vars number of variables
         * @param order polynomialmaximum degree
         * @param value polynomial constant term
         * @param monomial flag to detect the change of base (default false)
         */
        chebyshev_polynomial(const int &vars, const int &order, const T &value, const bool &monomial=false);

        ~chebyshev_polynomial();
        /******************************/
        /*ARITHMETIC OPERATIONS (+-*) */
        /******************************/
        /**
         * @brief operator + overload of operator sum.
         *
         * Implements the sum between two polynomials
         * @param other polynomial to be summed to the current polynomial
         * @return result of the sum operation
         */
         chebyshev_polynomial<T> operator+(const chebyshev_polynomial<T> &other) const;
         /**
          * @brief operator - overload of operator difference.
          *
          * Implements the dfference between two polynomials
          * @param other polynomial to be subtracted to the current polynomial
          * @return result of the difference operation
          */
         chebyshev_polynomial<T> operator-(const chebyshev_polynomial<T> &other) const;
         /**
          * @brief operator * overload of operator product.
          *
          * Implements the product between two polynomials. Stadard, DCT or monomial multiplication is performed according to the setting
          * @param other polynomial to be multiplied to the current polynomial
          * @return result of the multiplication
          */
        chebyshev_polynomial<T> operator*(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief operator / overload of divisin operator
         *
         * Implements the division between two polynomials.
         * @param other polynomial to be divided by
         * @return result of the division
         */
        chebyshev_polynomial<T> operator/(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief operator + sum of a constant term
         * @param other constant term
         * @return result of the sum
         */
        chebyshev_polynomial<T> operator+(const T& other) const;
        /**
         * @brief operator - difference with a constant term
         * @param other constant term
         * @return result of the difference
          */
        chebyshev_polynomial<T> operator-(const T& other) const;
        /**
         * @brief operator * multiplication by a constant term
         * @param other constant term
         * @return result of the multiplication
         */
        chebyshev_polynomial<T> operator*(const T& other) const;
        /**
         * @brief operator / division by a constant term
         * @param other constant term
         * @return result of the division
         */
        chebyshev_polynomial<T> operator/(const T& other) const;

        /******************************/
        /*UNARY OPERATORS             */
        /******************************/
        /**
         * @brief operator + return a polynomial with same coefficients
         * @return resulting polynomial
         */
        chebyshev_polynomial<T> operator+() const;
        /**
         * @brief operator - return a polynomial with opposite coefficients
         * @return resulting polynomial
         */
        chebyshev_polynomial<T> operator-() const;

        /******************************/
        /*ASSIGNEMENT (with operators)*/
        /******************************/
        /**
         * @brief operator = overload of assignment operator
         *
         * The coefficents of the polynomial on the rhs are copyed into the coefficients of the polynomial into the rhs
         * @param other polynomial rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator = overload of assignment operator
         *
         * The coefficients of the polynomial on the lhs are set to zero excluding the constant term that is initialized with
         * the rhs value
         * @param other value rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator=(const T &other);

        /**
         * @brief operator += overload of sum and assignment operator
         * @param other polynomial lhs
         * @return polynomial rhs
         */
        chebyshev_polynomial<T>& operator+=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator -= overload of differnece and assignment operator
         * @param other polynomial lhs
         * @return polynomial rhs
         */
        chebyshev_polynomial<T>& operator-=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator *= overload of product and assignment operator
         * @param other polynomial lhs
         * @return polynomial rhs
         */
        chebyshev_polynomial<T>& operator*=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator /= overload of division and assignment operator
         * @param other polynomial lhs
         * @return polynomial rhs
         */
        chebyshev_polynomial<T>& operator/=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator += overload of sum and assignment operator
         * @param other value rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator+=(const T& other);
        /**
         * @brief operator -= overload of difference and assignment operator
         * @param other value rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator-=(const T& other);
        /**
         * @brief operator *= overload of product and assignment operator
         * @param other value rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator*=(const T& other);
        /**
         * @brief operator /= overload of division and assignment operator
         * @param other value rhs
         * @return polynomial lhs
         */
        chebyshev_polynomial<T>& operator/=(const T& other);

        /******************************/
        /*COMPARISON                  */
        /******************************/
        /**
         * @brief operator == overload of comparison operator
         *
         * Two polynomials are equals when the set of coefficients are equals and are expressed in the same base
         * @param other polynomial rhs
         * @return boolean value
         */
        bool operator==(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief operator != overload of comparison operator
         *
         * Two polynomials are not equal if they differ for at least one coefficient
         * @param other polynomial rhs
         * @return boolean value
         */
        bool operator!=(const chebyshev_polynomial<T> &other) const;

        /**
         * @brief inv compute inverse of a polynomial by mean of chebychev approximation of 1/x and composition;
         * @param other polynomial to be inverted
         * @return inverse of the polynomial
         */
        chebyshev_polynomial<T> inv(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief direct_multiplication direct Chebyshev multipliation (highly inefficient)
         * @param x0 first factor of the multiplication
         * @param x1 second factor of the multiplication
         * @return resulting polynomial
         */
        chebyshev_polynomial<T> direct_multiplication(const chebyshev_polynomial<T> &x0, const chebyshev_polynomial<T> &x1) const;

    public:
        /******************************/
        /*EVALUATION & COMPOSITION    */
        /******************************/
        /**
         * @brief evaluate_1Dbase evaluate Chebyshev 1 dimensional base in a polynomial value and map it between [a,b]
         * @param other polynomial for evaluation
         * @param a lower interval bound
         * @param b upper interval bound
         * @return
         */
        static std::vector<chebyshev_polynomial<T> > evaluate_base1D(const chebyshev_polynomial<T> &other, const T &a, const T &b);

        /**
         * @brief composition compose the current polynomial with the vector of polynomial received as input
         * @param other vector of polynomials input for the composition
         */
        void composition(const std::vector<chebyshev_polynomial<T> > &other);

	/**
    * @brief evaluate_basis evaluate the multivariate chebyshev polynomial bases in a vector of constant values
    * @param x vector of constant values
    * @return the vector cntaining the evaluation of the basis in a lexicographic order
	*/
	std::vector<T> evaluate_basis(const std::vector<T> &x) const;
	/**
    * @brief evaluate evaluate the multivariate chebyshev polynomial in a vector of constant values
    * @param x vector of constant values
    * @return result of the polynomial evaluation
	*/
	T evaluate(const std::vector<T> &x) const;
	/**
     * @brief evaluate evaluate the univariate chebyshev polynomial in a constant value
     * @param x constant value
    * @return result of the polynomial evaluation
    */
    T evaluate(const T &x) const;

        /**
         * @brief map function for mapping the variables of the polynomial from [-1,1] to [a,b]
         * @param a vector defining the lower bound variables
         * @param b vector defining the upper bound variables
         */
    virtual void map(const std::vector<T> &a, const std::vector<T> &b);

    public:
        /******************************/
        /*CHANGE of BASIS             */
        /******************************/
    /**
     * @brief to_monomial_basis transformation to monomial base
     *
     * The method transforms the current polynomial into the monomial base.
     */
        void to_monomial_basis();
        /**
         * @brief from_monomial_basis transformation from monomial base
         *
         * The method transforms the current polynomial from the monomial base to its native one.
         */
        void from_monomial_basis();

        /**
         * @brief get_basis_name return a string that identifies the polynomial base
         *
         * It is used when polynomial is printed to std output
         * @return base identifier
         */
        std::string get_basis_name() const;

        /******************************/
        /*APPROXIMATION               */
        /******************************/
        /**
         * @brief approximation chebyshev approximation of a given function
         *
         * Approximate the univariate function f on the interval [a,b] up to specified degree
         * @param f function to be approximated
         * @param a lower interval bound
         * @param b upper interval bound
         * @param deg maximum degree of approximation
         * @return the vector of coefficients of the chebbyshev polynomial approximation of f
         */
        static std::vector<T> approximation(T (*f)(T x), const T &a, const T &b, const T &deg = chebyshev_polynomial<T>::MAX_DEGREE);

        /**
         * @brief approximation meta-function used for overloaded of elementary function.
         *
         * Approximate the function f and evaluates it in a polynomial
         * @param f function to be approximated
         * @param other polynomial for evaluation
         * @return the polynomial resulting from the evaluation of f in a polynomial
         */
        static chebyshev_polynomial<T> approximation(T (*f)(T x), const chebyshev_polynomial<T> &other);

    private:
        void initialize_t();
        std::vector<std::vector<int> > get_t() const {return m_t;}
        std::vector<std::vector<int> > m_t;
        //for fast 1D scalar evaluation
        T clenshaw(T x, int n) const;

	};
	
	template < class T>
	static chebyshev_polynomial<T> operator-(const T left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left)-right;
		return -right+left;
	}

	template < class T>
	static chebyshev_polynomial<T> operator-(const int left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)-right;
		return -right+((T) left);
	}

	template < class T>
	static chebyshev_polynomial<T> operator+(const T left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left)+right;
		return right+left;
	}

	template < class T>
	static chebyshev_polynomial<T> operator+(const int left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)+right;
		return right+((T) left);
	}

	template < class T>
	static chebyshev_polynomial<T> operator*(const T left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left)*right;
		return right*left;
	}

	template < class T>
	static chebyshev_polynomial<T> operator*(const int left, const chebyshev_polynomial<T> right){
		// return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)*right;
		return right*((T) left);
	}

	template < class T>
	static chebyshev_polynomial<T> operator/(const T left, const chebyshev_polynomial<T> right){
//        return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left, right.is_monomial_base())/right;
        return left*right.inv(right);
	}

	template < class T>
	static chebyshev_polynomial<T> operator/(const int left, const chebyshev_polynomial<T> right){
//		return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left, right.is_monomial_base())/right;
        return left*right.inv(right);
	}
}}



#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
