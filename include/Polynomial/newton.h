/******************************************************************************
 *                             NEWTON_POLYNOMIAL_H                            *
 *              Newton Polynomial Algebra of the SMART-UQ toolbox             *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#ifndef NEWTON_POLYNOMIAL_H_ 	
#define NEWTON_POLYNOMIAL_H_

#include "polynomial.h"

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
	class Newton_Polynomial: public Polynomial <T> {
	
	// for class template inheritance
	private:
		using Polynomial<T>::m_coeffs;
		using Polynomial<T>::m_degree;
		using Polynomial<T>::m_nvar;
		using Polynomial<T>::m_J;
		using Polynomial<T>::m_N;

	public:
		Newton_Polynomial(const int &vars, const int &order);
		//initialize newton polynomial equal to n1(x_i)
		Newton_Polynomial(const int &vars, const int &order, const int &i);
		//initialize a Newton polynomial equal to value
		Newton_Polynomial(const int &vars, const int &order, const T &value);

		static const int MAX_DEGREE = 100;
		// arithmetic operators
		Newton_Polynomial<T> operator+(const Newton_Polynomial<T> &other) const;
		Newton_Polynomial<T> operator-(const Newton_Polynomial<T> &other) const;
		Newton_Polynomial<T> operator*(const Newton_Polynomial<T> &other) const;
		Newton_Polynomial<T> operator/(const Newton_Polynomial<T> &other) const;
		Newton_Polynomial<T> operator+(const T& other) const;
		Newton_Polynomial<T> operator-(const T& other) const;
		Newton_Polynomial<T> operator*(const T& other) const;
		Newton_Polynomial<T> operator/(const T& other) const;
		// unary
		Newton_Polynomial<T> operator+() const;
		Newton_Polynomial<T> operator-() const;
		// assignment operator
		Newton_Polynomial<T>& operator=(const Newton_Polynomial<T> &other);
		Newton_Polynomial<T>& operator=(const T &other);
		// arithmetic operation and assignment
		Newton_Polynomial<T>& operator+=(const Newton_Polynomial<T> &other);
		Newton_Polynomial<T>& operator-=(const Newton_Polynomial<T> &other);
		Newton_Polynomial<T>& operator*=(const Newton_Polynomial<T> &other);
		Newton_Polynomial<T>& operator/=(const Newton_Polynomial<T> &other);
		Newton_Polynomial<T>& operator+=(const T& other);
		Newton_Polynomial<T>& operator-=(const T& other);
		Newton_Polynomial<T>& operator*=(const T& other);
		Newton_Polynomial<T>& operator/=(const T& other);

		//polynomial comparison (certainly) i.e. the comparison is true for every coefficient in the polynomial
		bool operator==(const Newton_Polynomial<T> &other) const;
		bool operator!=(const Newton_Polynomial<T> &other) const;

		// Newton_Polynomial<T> inv(const Newton_Polynomial<T> &other) const;
		// Newton_Polynomial<T> composition(const std::vector<Newton_Polynomial<T> > &other) const;

		// static Newton_Polynomial<T> direct_multiplication(const Newton_Polynomial<T> &x0, const Newton_Polynomial<T> &x1);

		// static std::vector<Newton_Polynomial<T> > evaluate_base(const Newton_Polynomial<T> &other, const T &a, const T &b);
		// static std::vector<T> cheb_approximation(T (*f)(T x), const T a, const T b);

		std::string get_basis_name() const;
		std::string get_name() const;

		T evaluate(const std::vector<T> &x) const;
		T evaluate(const T &x) const;

	public:
		std::vector<T> get_nodes() const {return m_nodes;}
		void set_nodes(std::vector<T> &nodes){
    		if(m_nodes.size()!=nodes.size()){
			std::cout<<"Provided inappropiate number of nodes."<<std::endl;
			exit(EXIT_FAILURE);
    		}
    	// 	for (int i=0;i<nodes.size();i++){
    	// 		if (fabs(nodes[i])>1.0){
					// std::cout<<"Nodes must belong to [-1,1]."<<std::endl;
					// exit(EXIT_FAILURE);
    	// 		}
    	// 	}
    		m_nodes=nodes;
		}
		
	private:
		void initialize_nodes();
		std::vector<T> m_nodes;

	private:
		T nested_multiplication(T x, int i) const;
	};
	
	template < class T>
	static Newton_Polynomial<T> operator-(const T left, const Newton_Polynomial<T> right){
		return -right+left;
	}

	template < class T>
	static Newton_Polynomial<T> operator-(const int left, const Newton_Polynomial<T> right){
		return -right+((T) left);
	}

	template < class T>
	static Newton_Polynomial<T> operator+(const T left, const Newton_Polynomial<T> right){
		return right+left;
	}

	template < class T>
	static Newton_Polynomial<T> operator+(const int left, const Newton_Polynomial<T> right){
		return right+((T) left);
	}

	template < class T>
	static Newton_Polynomial<T> operator*(const T left, const Newton_Polynomial<T> right){
		return right*left;
	}

	template < class T>
	static Newton_Polynomial<T> operator*(const int left, const Newton_Polynomial<T> right){
		return right*((T) left);
	}

	template < class T>
	static Newton_Polynomial<T> operator/(const T left, const Newton_Polynomial<T> right){
		return Newton_Polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
	}

	template < class T>
	static Newton_Polynomial<T> operator/(const int left, const Newton_Polynomial<T> right){
		return Newton_Polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)/right;
	}
}}



#endif /* NEWTON_POLYNOMIAL_H_ */
