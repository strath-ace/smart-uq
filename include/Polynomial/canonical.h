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

		canonical_polynomial(const int &vars, const int &order);
		//initialize a 1 degree univariate chebyshev polynomial of the corresponding variable [x1,x2,...]
		canonical_polynomial(const int &vars, const int &order, const int &i);
		//initialize a chebyshev polynomial with only the constant term
		canonical_polynomial(const int &vars, const int &order, const T &value);


		// arithmetic operators
		canonical_polynomial<T> operator+(const canonical_polynomial<T> &other) const;
		canonical_polynomial<T> operator-(const canonical_polynomial<T> &other) const;
		canonical_polynomial<T> operator*(const canonical_polynomial<T> &other) const;
		canonical_polynomial<T> operator/(const canonical_polynomial<T> &other) const;
		canonical_polynomial<T> operator+(const T& other) const;
		canonical_polynomial<T> operator-(const T& other) const;
		canonical_polynomial<T> operator*(const T& other) const;
		canonical_polynomial<T> operator/(const T& other) const;
		// unary
		canonical_polynomial<T> operator+() const;
		canonical_polynomial<T> operator-() const;
		// assignment operator
		canonical_polynomial<T>& operator=(const canonical_polynomial<T> &other);
		canonical_polynomial<T>& operator=(const T &other);
		// arithmetic operation and assignment
		canonical_polynomial<T>& operator+=(const canonical_polynomial<T> &other);
		canonical_polynomial<T>& operator-=(const canonical_polynomial<T> &other);
		canonical_polynomial<T>& operator*=(const canonical_polynomial<T> &other);
		canonical_polynomial<T>& operator/=(const canonical_polynomial<T> &other);
		canonical_polynomial<T>& operator+=(const T& other);
		canonical_polynomial<T>& operator-=(const T& other);
		canonical_polynomial<T>& operator*=(const T& other);
		canonical_polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		bool operator==(const canonical_polynomial<T> &other) const;
		bool operator!=(const canonical_polynomial<T> &other) const;

		canonical_polynomial<T> inv() const;
		canonical_polynomial<T> composition(const std::vector<canonical_polynomial<T> > &other) const;

		std::vector<canonical_polynomial<T> > evaluate_base(const T &a, const T &b) const;
		static std::vector<T> approximation_1d(T (*f)(T x), const T a, const T b, const int degree);
		void assign_from_chebyshev(const std::vector<T> cheb_coeffs);

		std::string get_basis_name() const;

		T evaluate(const std::vector<T> &x) const;
		T evaluate(const T &x) const;

		static void initialize_M(const int nvar, const int degree);
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
