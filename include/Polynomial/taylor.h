/******************************************************************************
 *                             taylor_polynomial_H                            *
 *               Taylor Polynomial Algebra of the SMART-UQ toolbox            *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
*/

#ifndef TAYLOR_H
#define TAYLOR_H

#include "base_polynomial.h"
#include "config.h"

using namespace std;

namespace smart{
namespace polynomial{
        template < class T >
    class taylor_polynomial: public base_polynomial <T> {

	// // for class template inheritance
	private:
	using base_polynomial<T>::m_name;
	using base_polynomial<T>::m_coeffs;
	using base_polynomial<T>::m_degree;
	using base_polynomial<T>::m_nvar;
	using base_polynomial<T>::m_J;
	using base_polynomial<T>::m_N;
	using base_polynomial<T>::m_manipulated_to_monomial;

        public:
        static const int MAX_DEGREE=100;
        /**
         * @brief taylor_polynomial
         * @param vars
         * @param order
         */
        taylor_polynomial(const int &vars, const int &order);
        /**
         * @brief taylor_polynomial
         * @param vars
         * @param order
         * @param i
         */
        taylor_polynomial(const int &vars, const int &order, const int &i);
        /**
         * @brief taylor_polynomial
         * @param vars
         * @param order
         * @param value
         */
        taylor_polynomial(const int &vars, const int &order, const T &value);
        /**
         * @brief taylor_polynomial
         * @param vars
         * @param order
         * @param i
         * @param a
         * @param b
         */
        taylor_polynomial(const int &vars, const int &order, const int &i, const T &a, const T &b);

        ~taylor_polynomial();
        /******************************/
        /*ARITHMETIC OPERATIONS (+-*) */
        /******************************/
        /**
         * @brief operator +
         * @param other
         * @return
         */
         taylor_polynomial<T> operator+(const taylor_polynomial<T> &other) const;
        /**
         * @brief operator -
         * @param other
         * @return
         */
         taylor_polynomial<T> operator-(const taylor_polynomial<T> &other) const;
        /**
         * @brief operator *
         * @param other
         * @return
         */
        taylor_polynomial<T> operator*(const taylor_polynomial<T> &other) const;
        /**
         * @brief operator /
         * @param other
         * @return
         */
        taylor_polynomial<T> operator/(const taylor_polynomial<T> &other) const;
        /**
         * @brief operator +
         * @param other
         * @return
         */
        taylor_polynomial<T> operator+(const T& other) const;
        /**
         * @brief operator -
         * @param other
         * @return
         */
        taylor_polynomial<T> operator-(const T& other) const;
        /**
         * @brief operator *
         * @param other
         * @return
         */
        taylor_polynomial<T> operator*(const T& other) const;
        /**
         * @brief operator /
         * @param other
         * @return
         */
        taylor_polynomial<T> operator/(const T& other) const;

        /**
         * @brief inv
         * @param other
         * @return
         */
        taylor_polynomial<T> inv(const taylor_polynomial<T> &other) const;

        /******************************/
        /*UNARY OPERATORS             */
        /******************************/
        /**
         * @brief operator +
         * @return
         */
        taylor_polynomial<T> operator+() const;
        /**
         * @brief operator -
         * @return
         */
        taylor_polynomial<T> operator-() const;

        /******************************/
        /*ASSIGNEMENT (with operators)*/
        /******************************/
        /**
         * @brief operator =
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator=(const taylor_polynomial<T> &other);
        /**
         * @brief operator =
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator=(const T &other);

        /**
         * @brief operator +=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator+=(const taylor_polynomial<T> &other);
        /**
         * @brief operator -=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator-=(const taylor_polynomial<T> &other);
        /**
         * @brief operator *=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator*=(const taylor_polynomial<T> &other);
        /**
         * @brief operator /=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator/=(const taylor_polynomial<T> &other);
        /**
         * @brief operator +=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator+=(const T& other);
        /**
         * @brief operator -=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator-=(const T& other);
        /**
         * @brief operator *=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator*=(const T& other);
        /**
         * @brief operator /=
         * @param other
         * @return
         */
        taylor_polynomial<T>& operator/=(const T& other);

        /******************************/
        /*COMPARISON                  */
        /******************************/
        /**
         * @brief operator ==
         * @param other
         * @return
         */
        bool operator==(const taylor_polynomial<T> &other) const;
        /**
         * @brief operator !=
         * @param other
         * @return
         */
        bool operator!=(const taylor_polynomial<T> &other) const;

    public:
        /******************************/
        /*EVALUATION & COMPOSITION    */
        /******************************/
        /**
         * @brief evaluate_1Dbase
         * @param other
         * @return
         */
        static std::vector<taylor_polynomial<T> > evaluate_base1D(const taylor_polynomial<T> &other);

        /**
         * @brief composition
         * @param other
         * @return
         */
        void composition(const std::vector<taylor_polynomial<T> > &other);

        /**
         * @brief evaluate_basis
         * @param x
         * @return
         */
        virtual std::vector<T> evaluate_basis(const std::vector<T> &x) const;
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

	virtual void map(const int &idx, const std::vector<T> &a, const std::vector<T> &b);

    public:
        /******************************/
        /*CHANGE of BASIS             */
        /******************************/
        /**
         * @brief to_monomial_basis
         * @return
         */
        void to_monomial_basis();
        /**
         * @brief from_monomial_basis
         * @param coeffs
         */
        void from_monomial_basis();

        /**
         * @brief get_basis_name
         * @return
         */
        std::string get_basis_name() const;

        /******************************/
        /*APPROXIMATION               */
        /******************************/
        /**
         * @brief approximation
         * @param a
         * @param b
         * @return
         */
        static std::vector<T> approximation(T (*f)(T x), const T &x0);

    private:
        T horner(T x, int i) const;

        };

	template < class T>
	static taylor_polynomial<T> operator-(const T left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), left)-right;
		return -right+left;
	}

	template < class T>
	static taylor_polynomial<T> operator-(const int left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)-right;
		return -right+((T) left);
	}

	template < class T>
	static taylor_polynomial<T> operator+(const T left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), left)+right;
		return right+left;
	}

	template < class T>
	static taylor_polynomial<T> operator+(const int left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)+right;
		return right+((T) left);
	}

	template < class T>
	static taylor_polynomial<T> operator*(const T left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), left)*right;
		return right*left;
	}

	template < class T>
	static taylor_polynomial<T> operator*(const int left, const taylor_polynomial<T> right){
		// return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)*right;
		return right*((T) left);
	}

	template < class T>
	static taylor_polynomial<T> operator/(const T left, const taylor_polynomial<T> right){
		return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), left)/right;
	}

	template < class T>
	static taylor_polynomial<T> operator/(const int left, const taylor_polynomial<T> right){
		return taylor_polynomial<T>(right.get_nvar(), right.get_degree(), (T) left)/right;
	}
}}


#endif // TAYLOR_H
