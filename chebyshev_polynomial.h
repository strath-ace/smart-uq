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
#include <numeric>
#include "utils.h"

using namespace std;

template < class T >
class Chebyshev_Polynomial{
public:
	Chebyshev_Polynomial(const int &vars, const int &order);
	//initialize a 1 degree univariate chebyshev polynomial of the corresponding variable [x1,x2,...]
	Chebyshev_Polynomial(const int &vars, const int &order, const int &i);
	//initialize a chebyshev polynomial with only the constant term
	Chebyshev_Polynomial(const int &vars, const int &order, const T &value);

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

	//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
	bool operator==(const Chebyshev_Polynomial<T> &other) const;
	bool operator!=(const Chebyshev_Polynomial<T> &other) const;

	friend ostream &operator<<(ostream &os, const Chebyshev_Polynomial<T> &poly) {
		std::vector<T> coeffs = poly.get_coeffs();
		int nvar = poly.get_nvar();
		int idx=0;

		os << "\t";
		for(int i=0; i<nvar; i++)
		    os << "x"<<i<<"\t";
		os << "\n";
		for(int deg=0; deg<=poly.get_degree(); deg++){
		    for(int i=0; i<poly.get_J()[poly.get_nvar()][deg]; i++){
			os << coeffs[idx] << "\t";
			std::vector<int> row = poly.get_row(i,deg);
			for(int j=0; j<row.size(); j++)
			    os<<row[j]<<"\t";
			os<<"\n";
			idx++;
		    }
		}
		return os;
	}

public:
	std::vector<T> get_coeffs() const {return m_coeffs;}
	void set_coeffs(std::vector<T> &coeffs){m_coeffs=coeffs;}
	int get_degree() const {return m_degree;}
	int get_nvar() const {return m_nvar;}
	std::vector<std::vector<int> > get_J() const {return m_J;}

	int get_ncoeffs() const {return m_coeffs.size();}

private:
	//BEGIN A&V
	Chebyshev_Polynomial(const T& l, const T& r,const bool& b);
	Chebyshev_Polynomial<T> add(const Chebyshev_Polynomial<T> &A,const Chebyshev_Polynomial<T> &B) const;
	Chebyshev_Polynomial<T> sub(const Chebyshev_Polynomial<T> &A,const Chebyshev_Polynomial<T> &B) const;
	Chebyshev_Polynomial<T> mult(const Chebyshev_Polynomial<T> &A, const Chebyshev_Polynomial<T> &B) const;
	//END

private:
	//polynomial representation variables
	void initialize_J();
	void initialize_N();
	void initialize_t();
	std::vector<int> get_row(const int &idx, const int &deg) const;
	int get_idx(const std::vector<int> &k) const;

private:
	std::vector<T> m_coeffs;
	int m_degree;
	int m_nvar;
	std::vector<std::vector<int> > m_J, m_N, m_t;
};

template < class T>
static Chebyshev_Polynomial<T> operator-(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)-right;
}

template < class T>
static Chebyshev_Polynomial<T> operator-(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)-right;
}

template < class T>
static Chebyshev_Polynomial<T> operator+(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)+right;
}

template < class T>
static Chebyshev_Polynomial<T> operator+(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)+right;
}

template < class T>
static Chebyshev_Polynomial<T> operator*(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)*right;
}

template < class T>
static Chebyshev_Polynomial<T> operator*(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)*right;
}

template < class T>
static Chebyshev_Polynomial<T> operator/(const T left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
}

template < class T>
static Chebyshev_Polynomial<T> operator/(const int left, const Chebyshev_Polynomial<T> right){
	return Chebyshev_Polynomial<T>(right.get_nvar(), right.get_degree(), left)/-right;
}


#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
