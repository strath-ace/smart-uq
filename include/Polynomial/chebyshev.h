/******************************************************************************
 *                           chebyshev_polynomial_H                           *
 *            Chebyshev Polynomial Algebra of the SMART-UQ toolbox            *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include "polynomial.h"
#include "config.h"

#ifdef CHEBYSHEV_DCT_MULTIPLICATION
#include <fftw3.h>
#endif

using namespace std;

namespace smart{
namespace polynomial_algebra{
	template < class T >
    class chebyshev_polynomial: public polynomial <T> {
	
	// // for class template inheritance
	private:
    using polynomial<T>::m_name;
    using polynomial<T>::m_coeffs;
    using polynomial<T>::m_degree;
    using polynomial<T>::m_nvar;
    using polynomial<T>::m_J;
    using polynomial<T>::m_N;

	public:
        static const int MAX_DEGREE=100;
		/**
		 * @brief chebyshev_polynomial
		 * @param vars
		 * @param order
		 */
		chebyshev_polynomial(const int &vars, const int &order);
		/**
		 * @brief chebyshev_polynomial
		 * @param vars
		 * @param order
		 * @param i
		 */
		chebyshev_polynomial(const int &vars, const int &order, const int &i);
		/**
		 * @brief chebyshev_polynomial
		 * @param vars
		 * @param order
		 * @param value
		 */
		chebyshev_polynomial(const int &vars, const int &order, const T &value);

        /******************************/
        /*ARITHMETIC OPERATIONS (+-*) */
        /******************************/
        /**
         * @brief operator +
         * @param other
         * @return
         */
         chebyshev_polynomial<T> operator+(const chebyshev_polynomial<T> &other) const{return polynomial<T>::operator +(other);}
        /**
         * @brief operator -
         * @param other
         * @return
         */
         chebyshev_polynomial<T> operator-(const chebyshev_polynomial<T> &other) const{return polynomial<T>::operator -(other);}
        /**
         * @brief operator *
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator*(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief operator /
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator/(const chebyshev_polynomial<T> &other) const;
        /**
         * @brief operator +
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator+(const T& other) const{return polynomial<T>::operator +(other);}
        /**
         * @brief operator -
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator-(const T& other) const{return polynomial<T>::operator -(other);}
        /**
         * @brief operator *
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator*(const T& other) const{return polynomial<T>::operator *(other);}
        /**
         * @brief operator /
         * @param other
         * @return
         */
        chebyshev_polynomial<T> operator/(const T& other) const{return polynomial<T>::operator /(other);}

        /******************************/
        /*UNARY OPERATORS             */
        /******************************/
        /**
         * @brief operator +
         * @return
         */
        chebyshev_polynomial<T> operator+() const{return polynomial<T>::operator +();}
        /**
         * @brief operator -
         * @return
         */
        chebyshev_polynomial<T> operator-() const{return polynomial<T>::operator -();}

        /******************************/
        /*ASSIGNEMENT (with operators)*/
        /******************************/
        /**
         * @brief operator =
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator=(const chebyshev_polynomial<T> &other){polynomial<T>::operator =(other); return *this;}
        /**
         * @brief operator =
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator=(const T &other){polynomial<T>::operator =(other); return *this;}

        /**
         * @brief operator +=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator+=(const chebyshev_polynomial<T> &other){polynomial<T>::operator +=(other); return *this;}
        /**
         * @brief operator -=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator-=(const chebyshev_polynomial<T> &other){polynomial<T>::operator -=(other); return *this;}
        /**
         * @brief operator *=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator*=(const chebyshev_polynomial<T> &other){polynomial<T>::operator *=(other); return *this;}
        /**
         * @brief operator /=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator/=(const chebyshev_polynomial<T> &other);
        /**
         * @brief operator +=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator+=(const T& other){polynomial<T>::operator +=(other); return *this;}
        /**
         * @brief operator -=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator-=(const T& other){polynomial<T>::operator -=(other); return *this;}
        /**
         * @brief operator *=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator*=(const T& other){polynomial<T>::operator *=(other); return *this;}
        /**
         * @brief operator /=
         * @param other
         * @return
         */
        chebyshev_polynomial<T>& operator/=(const T& other){polynomial<T>::operator /=(other); return *this;}

        /******************************/
        /*COMPARISON                  */
        /******************************/
        /**
         * @brief operator ==
         * @param other
         * @return
         */
        bool operator==(const chebyshev_polynomial<T> &other) const{return polynomial<T>::operator ==(other);}
        /**
         * @brief operator !=
         * @param other
         * @return
         */
        bool operator!=(const chebyshev_polynomial<T> &other) const{return polynomial<T>::operator !=(other);}

        /******************************/
        /*I/O OPERATOR                */
        /******************************/
        /**
         * @brief operator <<
         * @param os
         * @param poly
         * @return
         */
        friend ostream &operator<<(ostream &os, const chebyshev_polynomial<T> &poly){return polynomial<T>::operator <<(os,poly);}

        /**
		 * @brief inv
		 * @param other
		 * @return
		 */
		chebyshev_polynomial<T> inv(const chebyshev_polynomial<T> &other) const;

		/**
		 * @brief direct_multiplication
		 * @param x0
		 * @param x1
		 * @return
		 */
        chebyshev_polynomial<T> direct_multiplication(const chebyshev_polynomial<T> &x0, const chebyshev_polynomial<T> &x1) const;

    public:
        /******************************/
        /*EVALUATION & COMPOSITION    */
        /******************************/
        /**
         * @brief evaluate_1Dbase
         * @param other
         * @param a
         * @param b
         * @return
         */
        std::vector<chebyshev_polynomial<T> > evaluate_base1D(const chebyshev_polynomial<T> &other, const T &a, const T &b) const;
        /**
         * @brief evaluate_1Dbase
         * @param other
         * @return
         */
        std::vector<T> evaluate_base1D(const T &other, const T &a, const T &b) const;
        /**
         * @brief composition
         * @param other
         * @return
         */
        chebyshev_polynomial<T> composition(const std::vector<chebyshev_polynomial<T> > &other) const;

		/**
		 * @brief evaluate
		 * @param x
		 * @return
		 */
		T evaluate(const std::vector<T> &x) const;

    public:
        /******************************/
        /*CHANGE of BASIS             */
        /******************************/
        /**
         * @brief to_monomial_basis
         * @return
         */
        std::vector<T> to_monomial_basis() const;
        /**
         * @brief to_monomial_basis
         * @param coeffs
         * @return
         */
        std::vector<T> to_monomial_basis(const std::vector<T> &coeffs) const;
        /**
         * @brief from_monomial_basis
         * @param coeffs
         */
        void from_monomial_basis(const std::vector<T> &coeffs) const;

        /**
         * @brief get_basis_name
         * @return
         */
        std::string get_basis_name() const;

        /******************************/
        /*APPROXIMATION               */
        /******************************/
        /**
         * @brief cheb_approximation
         * @param a
         * @param b
         * @return
         */
        static std::vector<T> approximation(T (*f)(T x), const T &a, const T &b);

        /**
         * @brief cheb_approximation
         * @param a
         * @param b
         * @param other
         * @return
         */
        static chebyshev_polynomial<T> approximation(T (*f)(T x), const chebyshev_polynomial<T> &other);

    private:
        void initialize_t();
		std::vector<std::vector<int> > get_t() const {return m_t;}
		std::vector<std::vector<int> > m_t;



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
		return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
	}

	template < class T>
	static chebyshev_polynomial<T> operator/(const int left, const chebyshev_polynomial<T> right){
		return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)/right;
	}
}}



#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
