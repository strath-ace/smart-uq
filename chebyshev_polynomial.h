/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "utils.h"

using namespace std;

template < class T >
class Chebyshev_Polynomial{
public:
	//empty interval
	Chebyshev_Polynomial(int vars, int order);

	//deconstructor
	~Chebyshev_Polynomial(){}

	//arithmetic operators
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

	//interval comparison (certainly) i.e. the comparison is true for every pair of element of the two intervals
	bool operator==(const Chebyshev_Polynomial<T> &other) const;
	bool operator!=(const Chebyshev_Polynomial<T> &other) const;
	bool operator<(const Chebyshev_Polynomial<T> &other) const;
	bool operator<=(const Chebyshev_Polynomial<T> &other) const;
	bool operator>(const Chebyshev_Polynomial<T> &other) const;
	bool operator>=(const Chebyshev_Polynomial<T> &other) const;

	friend ostream &operator<<(ostream &os, const Chebyshev_Polynomial<T> &i) {
		std::vector<T> coeffs = i.get_coeffs();
		int n = i.get_ncoeffs();
		int nvar = i.get_nvar();

		os << "\t ";
		for(int i=0; i<nvar; i++)
		    os << "x"<<i<<"\t";
		os << "\n";
		for(int i=0; i<n; i++){
		    os<<coeffs[i]<<"\n";
		}
		return os;
	}

public:
	std::vector<T>& get_coeffs(){return m_coeffs;}
	void set_coeffs(std::vector<T> &coeffs){m_coeffs=coeffs;}
	int get_degree(){return m_degree;}
	int get_nvar(){return m_nvar;}

	int get_ncoeffs(){return m_coeffs.size();}

private:
	//BEGIN A&V
	Chebyshev_Polynomial(const T& l, const T& r,const bool& b);
	Chebyshev_Polynomial<T> add(const Chebyshev_Polynomial<T> &A,const Chebyshev_Polynomial<T> &B) const;
	Chebyshev_Polynomial<T> sub(const Chebyshev_Polynomial<T> &A,const Chebyshev_Polynomial<T> &B) const;
	Chebyshev_Polynomial<T> mult(const Chebyshev_Polynomial<T> &A, const Chebyshev_Polynomial<T> &B) const;
	//END

private:
	std::vector<T> m_coeffs;
	const int m_degree;
	const int m_nvar;
};

template < class T>
static Chebyshev_Polynomial<T> operator-(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)-right;
}

template < class T>
static Chebyshev_Polynomial<T> operator-(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)-right;
}

template < class T>
static Chebyshev_Polynomial<T> operator+(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)+right;
}

template < class T>
static Chebyshev_Polynomial<T> operator+(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)+right;
}

template < class T>
static Chebyshev_Polynomial<T> operator*(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)*right;
}

template < class T>
static Chebyshev_Polynomial<T> operator*(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)*right;
}

template < class T>
static Chebyshev_Polynomial<T> operator/(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)/right;
}

template < class T>
static Chebyshev_Polynomial<T> operator/(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(left)/-right;
}


#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
