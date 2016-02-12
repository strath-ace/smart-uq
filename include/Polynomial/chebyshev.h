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

#include "base_polynomial.h"
#include "config.h"

#ifdef CHEBYSHEV_DCT_MULTIPLICATION
#include <fftw3.h>
#endif

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
	class chebyshev_polynomial: public base_polynomial <T> {
	
	// // for class template inheritance
	private:
	using base_polynomial<T>::m_name;
	using base_polynomial<T>::m_coeffs;
	using base_polynomial<T>::m_degree;
	using base_polynomial<T>::m_nvar;
	using base_polynomial<T>::m_J;
	using base_polynomial<T>::m_N;


	public:
		static const int MAX_DEGREE = 100;

		chebyshev_polynomial(const int &vars, const int &order);
		//initialize a 1 degree univariate chebyshev polynomial of the corresponding variable [x1,x2,...]
		chebyshev_polynomial(const int &vars, const int &order, const int &i);
		//initialize a chebyshev polynomial with only the constant term
		chebyshev_polynomial(const int &vars, const int &order, const T &value);

		// arithmetic operators
		chebyshev_polynomial<T> operator+(const chebyshev_polynomial<T> &other) const;
		chebyshev_polynomial<T> operator-(const chebyshev_polynomial<T> &other) const;
		chebyshev_polynomial<T> operator*(const chebyshev_polynomial<T> &other) const;
		chebyshev_polynomial<T> operator/(const chebyshev_polynomial<T> &other) const;
		chebyshev_polynomial<T> operator+(const T& other) const;
		chebyshev_polynomial<T> operator-(const T& other) const;
		chebyshev_polynomial<T> operator*(const T& other) const;
		chebyshev_polynomial<T> operator/(const T& other) const;
		//unary
		chebyshev_polynomial<T> operator+() const;
		chebyshev_polynomial<T> operator-() const;
		//assignment operator
		chebyshev_polynomial<T>& operator=(const chebyshev_polynomial<T> &other);
		chebyshev_polynomial<T>& operator=(const T &other);
		//arithmetic operation and assignment
		chebyshev_polynomial<T>& operator+=(const chebyshev_polynomial<T> &other);
		chebyshev_polynomial<T>& operator-=(const chebyshev_polynomial<T> &other);
		chebyshev_polynomial<T>& operator*=(const chebyshev_polynomial<T> &other);
		chebyshev_polynomial<T>& operator/=(const chebyshev_polynomial<T> &other);
		chebyshev_polynomial<T>& operator+=(const T& other);
		chebyshev_polynomial<T>& operator-=(const T& other);
		chebyshev_polynomial<T>& operator*=(const T& other);
		chebyshev_polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		bool operator==(const chebyshev_polynomial<T> &other) const;
		bool operator!=(const chebyshev_polynomial<T> &other) const;

		chebyshev_polynomial<T> inv(const chebyshev_polynomial<T> &other) const;
		chebyshev_polynomial<T> composition(const std::vector<chebyshev_polynomial<T> > &other) const;

		static chebyshev_polynomial<T> direct_multiplication(const chebyshev_polynomial<T> &x0, const chebyshev_polynomial<T> &x1);

		static std::vector<chebyshev_polynomial<T> > evaluate_base(const chebyshev_polynomial<T> &other, const T &a, const T &b);
		static std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b);

		std::string get_basis_name() const;

		T evaluate(const std::vector<T> &x) const;
		T evaluate(const T &x) const;
	
	private:
		void initialize_t();
		std::vector<std::vector<int> > get_t() const {return m_t;}
		std::vector<std::vector<int> > m_t;

	private:
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
		return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
	}

	template < class T>
	static chebyshev_polynomial<T> operator/(const int left, const chebyshev_polynomial<T> right){
		return chebyshev_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)/right;
	}
}}



#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
