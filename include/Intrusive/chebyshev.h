/******************************************************************************
 *                           CHEBYSHEV_POLYNOMIAL_H                           *
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
#include <fftw3.h>

using namespace std;

namespace smart{
namespace intrusive{
	template < class T >
	class Chebyshev_Polynomial: public Polynomial <T> {
	public:
		static const int MAX_DEGREE = 100;

		using Polynomial<T>::Polynomial;

		// arithmetic operators
		Chebyshev_Polynomial<T> operator+(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> operator-(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> operator*(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> operator/(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> operator+(const T& other) const;
		Chebyshev_Polynomial<T> operator-(const T& other) const;
		Chebyshev_Polynomial<T> operator*(const T& other) const;
		Chebyshev_Polynomial<T> operator/(const T& other) const;
		//assignment operator
		Chebyshev_Polynomial<T>& operator=(const Chebyshev_Polynomial<T> &other);
		Chebyshev_Polynomial<T>& operator=(const T &other);
		//arithmetic operation and assignment
		Chebyshev_Polynomial<T>& operator+=(const Chebyshev_Polynomial<T> &other);
		Chebyshev_Polynomial<T>& operator-=(const Chebyshev_Polynomial<T> &other);
		Chebyshev_Polynomial<T>& operator*=(const Chebyshev_Polynomial<T> &other);
		Chebyshev_Polynomial<T>& operator/=(const Chebyshev_Polynomial<T> &other);
		Chebyshev_Polynomial<T>& operator+=(const T& other);
		Chebyshev_Polynomial<T>& operator-=(const T& other);
		Chebyshev_Polynomial<T>& operator*=(const T& other);
		Chebyshev_Polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		bool operator==(const Chebyshev_Polynomial<T> &other) const;
		bool operator!=(const Chebyshev_Polynomial<T> &other) const;

		Chebyshev_Polynomial<T> inv(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> composition(const std::vector<Chebyshev_Polynomial<T> > &other) const;
		std::string get_basis() const;

		static std::vector<Chebyshev_Polynomial<T> > evaluate_base(const Chebyshev_Polynomial<T> &other, const T &a, const T &b);
		static std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b);

	public:
	using Polynomial<T>::m_coeffs;
	using Polynomial<T>::m_degree;
	using Polynomial<T>::m_nvar;

	private:
	using Polynomial<T>::m_J;
	using Polynomial<T>::m_N;
	using Polynomial<T>::m_t;
	};
	
	template < class T>
	static Chebyshev_Polynomial<T> operator-(const T left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)-right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator-(const int left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)-right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator+(const T left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)+right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator+(const int left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)+right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator*(const T left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)*right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator*(const int left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)*right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator/(const T left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
	}

	template < class T>
	static Chebyshev_Polynomial<T> operator/(const int left, const Chebyshev_Polynomial<T> right){
		return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)/right;
	}
}}



#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
