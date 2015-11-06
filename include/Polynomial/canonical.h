/******************************************************************************
 *                             Canonical_POLYNOMIAL_H                            *
 *              Canonical Polynomial Algebra of the SMART-UQ toolbox             *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#ifndef CANONICAL_POLYNOMIAL_H_ 	
#define CANONICAL_POLYNOMIAL_H_

#include "polynomial.h"

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
	class Canonical_Polynomial: public Polynomial <T> {
	
	// for class template inheritance
	public:
		using Polynomial<T>::Polynomial;
	// for class template inheritance
	private:
		using Polynomial<T>::m_coeffs;
		using Polynomial<T>::m_degree;
		using Polynomial<T>::m_nvar;
		using Polynomial<T>::m_J;
		using Polynomial<T>::m_N;

	public:

		static const int MAX_DEGREE = 100;

		// arithmetic operators
		Canonical_Polynomial<T> operator+(const Canonical_Polynomial<T> &other) const;
		Canonical_Polynomial<T> operator-(const Canonical_Polynomial<T> &other) const;
		Canonical_Polynomial<T> operator*(const Canonical_Polynomial<T> &other) const;
		Canonical_Polynomial<T> operator/(const Canonical_Polynomial<T> &other) const;
		Canonical_Polynomial<T> operator+(const T& other) const;
		Canonical_Polynomial<T> operator-(const T& other) const;
		Canonical_Polynomial<T> operator*(const T& other) const;
		Canonical_Polynomial<T> operator/(const T& other) const;
		// unary
		Canonical_Polynomial<T> operator+() const;
		Canonical_Polynomial<T> operator-() const;
		// assignment operator
		Canonical_Polynomial<T>& operator=(const Canonical_Polynomial<T> &other);
		Canonical_Polynomial<T>& operator=(const T &other);
		// arithmetic operation and assignment
		Canonical_Polynomial<T>& operator+=(const Canonical_Polynomial<T> &other);
		Canonical_Polynomial<T>& operator-=(const Canonical_Polynomial<T> &other);
		Canonical_Polynomial<T>& operator*=(const Canonical_Polynomial<T> &other);
		Canonical_Polynomial<T>& operator/=(const Canonical_Polynomial<T> &other);
		Canonical_Polynomial<T>& operator+=(const T& other);
		Canonical_Polynomial<T>& operator-=(const T& other);
		Canonical_Polynomial<T>& operator*=(const T& other);
		Canonical_Polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		bool operator==(const Canonical_Polynomial<T> &other) const;
		bool operator!=(const Canonical_Polynomial<T> &other) const;

		Canonical_Polynomial<T> inv() const;
		Canonical_Polynomial<T> composition(const std::vector<Canonical_Polynomial<T> > &other) const;

		// static Canonical_Polynomial<T> direct_multiplication(const Canonical_Polynomial<T> &x0, const Canonical_Polynomial<T> &x1);

		std::vector<Canonical_Polynomial<T> > evaluate_base(const T &a, const T &b) const;
		static std::vector<T> approximation_1d(T (*f)(T x), const T a, const T b, const int degree);

		std::string get_basis_name() const;
		std::string get_name() const;

		T evaluate(const std::vector<T> &x) const;
		T evaluate(const T &x) const;

	private:
		T horner(T x, int i) const;

	};
	
	template < class T>
	static Canonical_Polynomial<T> operator-(const T left, const Canonical_Polynomial<T> right){
		return -right+left;
	}

	template < class T>
	static Canonical_Polynomial<T> operator-(const int left, const Canonical_Polynomial<T> right){
		return -right+((T) left);
	}

	template < class T>
	static Canonical_Polynomial<T> operator+(const T left, const Canonical_Polynomial<T> right){
		return right+left;
	}

	template < class T>
	static Canonical_Polynomial<T> operator+(const int left, const Canonical_Polynomial<T> right){
		return right+((T) left);
	}

	template < class T>
	static Canonical_Polynomial<T> operator*(const T left, const Canonical_Polynomial<T> right){
		return right*left;
	}

	template < class T>
	static Canonical_Polynomial<T> operator*(const int left, const Canonical_Polynomial<T> right){
		return right*((T) left);
	}

	template < class T>
	static Canonical_Polynomial<T> operator/(const T left, const Canonical_Polynomial<T> right){
		return left*right.inv();
	}

	template < class T>
	static Canonical_Polynomial<T> operator/(const int left, const Canonical_Polynomial<T> right){
		return ((T) left)*right.inv();
	}
}}



#endif /* CANONICAL_POLYNOMIAL_H_ */
