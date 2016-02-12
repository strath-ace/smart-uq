/******************************************************************************
 *                             canonical_polynomial_H                            *
 *              Canonical Polynomial Algebra of the SMART-UQ toolbox             *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#ifndef CANONICAL_POLYNOMIAL_H_
#define CANONICAL_POLYNOMIAL_H_

#include "base_polynomial.h"

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
	class canonical_polynomial: public base_polynomial <T> {

	// for class template inheritance
	private:
		using base_polynomial<T>::m_coeffs;
		using base_polynomial<T>::m_degree;
		using base_polynomial<T>::m_nvar;
		using base_polynomial<T>::m_J;
		using base_polynomial<T>::m_N;
		using base_polynomial<T>::m_name;

	public:

		static const int MAX_DEGREE = 100;

		/**
		 * @brief canonical_polynomial
		 * @param vars
		 * @param order
		 */
		canonical_polynomial(const int &vars, const int &order);
		/**
		 * @brief canonical_polynomial
		 * @param vars
		 * @param order
		 * @param i
		 */
		canonical_polynomial(const int &vars, const int &order, const int &i);
		/**
		 * @brief canonical_polynomial
		 * @param vars
		 * @param order
		 * @param value
		 */
		canonical_polynomial(const int &vars, const int &order, const T &value);


		// arithmetic operators
		/**
		 * @brief operator +
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator+(const canonical_polynomial<T> &other) const;
		/**
		 * @brief operator -
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator-(const canonical_polynomial<T> &other) const;
		/**
		 * @brief operator *
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator*(const canonical_polynomial<T> &other) const;
		/**
		 * @brief operator /
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator/(const canonical_polynomial<T> &other) const;
		/**
		 * @brief operator +
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator+(const T& other) const;
		/**
		 * @brief operator -
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator-(const T& other) const;
		/**
		 * @brief operator *
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator*(const T& other) const;
		/**
		 * @brief operator /
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> operator/(const T& other) const;

		// unary
		/**
		 * @brief operator +
		 * @return
		 */
		canonical_polynomial<T> operator+() const;
		/**
		 * @brief operator -
		 * @return
		 */
		canonical_polynomial<T> operator-() const;

		// assignment operator
		/**
		 * @brief operator =
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator=(const canonical_polynomial<T> &other);
		/**
		 * @brief operator =
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator=(const T &other);
		// arithmetic operation and assignment
		/**
		 * @brief operator +=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator+=(const canonical_polynomial<T> &other);
		/**
		 * @brief operator -=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator-=(const canonical_polynomial<T> &other);
		/**
		 * @brief operator *=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator*=(const canonical_polynomial<T> &other);
		/**
		 * @brief operator /=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator/=(const canonical_polynomial<T> &other);
		/**
		 * @brief operator +=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator+=(const T& other);
		/**
		 * @brief operator -=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator-=(const T& other);
		/**
		 * @brief operator *=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator*=(const T& other);
		/**
		 * @brief operator /=
		 * @param other
		 * @return
		 */
		canonical_polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		/**
		 * @brief operator ==
		 * @param other
		 * @return
		 */
		bool operator==(const canonical_polynomial<T> &other) const;
		/**
		 * @brief operator !=
		 * @param other
		 * @return
		 */
		bool operator!=(const canonical_polynomial<T> &other) const;

		/**
		 * @brief inv
		 * @return
		 */
		canonical_polynomial<T> inv() const;
		/**
		 * @brief composition
		 * @param other
		 * @return
		 */
		canonical_polynomial<T> composition(const std::vector<canonical_polynomial<T> > &other) const;

		/**
		 * @brief evaluate_base
		 * @param a
		 * @param b
		 * @return
		 */
		std::vector<canonical_polynomial<T> > evaluate_base(const T &a, const T &b) const;
		/**
		 * @brief approximation_1d
		 * @param a
		 * @param b
		 * @param degree
		 * @return
		 */
		static std::vector<T> approximation_1d(T (*f)(T x), const T a, const T b, const int degree);
		/**
		 * @brief assign_from_chebyshev
		 * @param cheb_coeffs
		 */
		void assign_from_chebyshev(const std::vector<T> cheb_coeffs);

		/**
		 * @brief get_basis_name
		 * @return
		 */
		std::string get_basis_name() const;

		/**
		 * @brief evaluate
		 * @param x
		 * @return
		 */
		T evaluate(const std::vector<T> &x) const;
		/**
		 * @brief evaluate
		 * @param x
		 * @return
		 */
		T evaluate(const T &x) const;

		/**
		 * @brief initialize_M
		 * @param nvar
		 * @param degree
		 */
		static void initialize_M(const int nvar, const int degree);
		/**
		 * @brief delete_M
		 */
		static void delete_M();

	private:
		T horner(T x, int i) const;
		static std::vector<int> m_M;
		static int m_Mnvar, m_Mdegree;
	};
	
	template < class T>
	static canonical_polynomial<T> operator-(const T left, const canonical_polynomial<T> right){
		return -right+left;
	}

	template < class T>
	static canonical_polynomial<T> operator-(const int left, const canonical_polynomial<T> right){
		return -right+((T) left);
	}

	template < class T>
	static canonical_polynomial<T> operator+(const T left, const canonical_polynomial<T> right){
		return right+left;
	}

	template < class T>
	static canonical_polynomial<T> operator+(const int left, const canonical_polynomial<T> right){
		return right+((T) left);
	}

	template < class T>
	static canonical_polynomial<T> operator*(const T left, const canonical_polynomial<T> right){
		return right*left;
	}

	template < class T>
	static canonical_polynomial<T> operator*(const int left, const canonical_polynomial<T> right){
		return right*((T) left);
	}

	template < class T>
	static canonical_polynomial<T> operator/(const T left, const canonical_polynomial<T> right){
		return left*right.inv();
	}

	template < class T>
	static canonical_polynomial<T> operator/(const int left, const canonical_polynomial<T> right){
		return ((T) left)*right.inv();
	}
}}



#endif /* CANONICAL_POLYNOMIAL_H_ */
