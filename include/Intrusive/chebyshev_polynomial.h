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

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <iomanip>
#include "utils.h"
#include <fftw3.h>

using namespace std;

namespace smart{
namespace intrusive{
	template < class T >
	class Chebyshev_Polynomial{
	public:
		static const int MAX_DEGREE = 100;

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

		Chebyshev_Polynomial<T> inv(const Chebyshev_Polynomial<T> &other) const;
		Chebyshev_Polynomial<T> composition(const std::vector<Chebyshev_Polynomial<T> > &other) const;

		static std::vector<Chebyshev_Polynomial<T> > evaluate_base(const Chebyshev_Polynomial<T> &other, const T &a, const T &b);
		static std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b);

		friend ostream &operator<<(ostream &os, const Chebyshev_Polynomial<T> &poly) {
			std::vector<T> coeffs = poly.get_coeffs();
			int nvar = poly.get_nvar();
			int idx=0;

			os <<std::setfill(' ')<<setw(16);
			for(int i=0; i<nvar; i++)
		    		os << "T(x"<<i<<")\t";
			os << "\n";
			for(int deg=0; deg<=poly.get_degree(); deg++){
	    		for(int i=0; i<poly.get_J()[poly.get_nvar()][deg]; i++){
					os <<left<<setw(16)<<coeffs[idx];
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
       		//getter and setters
		std::vector<T> get_coeffs() const {return m_coeffs;}
		void set_coeffs(std::vector<T> &coeffs){
	    		if(m_coeffs.size()!=coeffs.size()){
				std::cout<<"Coefficients vectors don't have the same lenght"<<std::endl;
				exit(EXIT_FAILURE);
	    		}
	    		m_coeffs=coeffs;
		}
		void set_coeffs(const int &idx, const T &value){m_coeffs[idx]=value;}
		int get_degree() const {return m_degree;}
		int get_nvar() const {return m_nvar;}
		std::vector<std::vector<int> > get_J() const {return m_J;}
		std::vector<std::vector<int> > get_N() const {return m_N;}
		std::vector<int> get_row(const int &idx, const int &deg) const;
		std::vector<std::vector<int> > get_t() const {return m_t;}
		int get_idx(const std::vector<int> &k) const;

		std::vector<T> get_range() const{
	    		std::vector<T> range(2,0);
	    		T constant = m_coeffs[0];
	    		for(int i=1; i<m_coeffs.size(); i++)
				range[1] += fabs(m_coeffs[i]);
	    		range[0] = -range[1]+constant;
	    		range[1] += constant;
	    		return range;
		}

	private:
		//polynomial representation variables
		void initialize_J();
		void initialize_N();
		void initialize_t();


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
}
}



#endif /* CHEBYSHEV_POLYNOMIAL_H_ */
