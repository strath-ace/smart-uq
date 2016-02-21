/******************************************************************************
 *                                  polynomial_H                              *
 *              Polynomial Algebra superclass of the SMART-UQ toolbox         *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde-----------------
------------------ e-mail: carlos.ortega@strath.ac.uk ------------------------
--------------------------- Author: Carlos Ortega ----------------------------
*/

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <numeric>
#include <iomanip>
#include "exception.h"
#include "constants.h"
#include "utils.h"

using namespace std;

namespace smart{
namespace polynomial_algebra{
	template < class T >
    class polynomial{

	public:
        /**************/
        /*CONSTRUCTORS*/
        /**************/
		/**
         * @brief polynomial
		 * @param vars
		 * @param order
		 */
        polynomial(const int &vars, const int &order);
		/**
         * @brief polynomial
		 * @param vars
		 * @param order
		 * @param i
		 */
        polynomial(const int &vars, const int &order, const int &i);
		/**
         * @brief polynomial
		 * @param vars
		 * @param order
		 * @param value
		 */
        polynomial(const int &vars, const int &order, const T &value);


        /******************************/
        /*ARITHMETIC OPERATIONS (+-*) */
        /******************************/
        /**
		 * @brief operator +
		 * @param other
		 * @return
		 */
        polynomial<T> operator+(const polynomial<T> &other) const;
		/**
		 * @brief operator -
		 * @param other
		 * @return
		 */
        polynomial<T> operator-(const polynomial<T> &other) const;
		/**
		 * @brief operator *
		 * @param other
		 * @return
		 */
        polynomial<T> operator*(const polynomial<T> &other) const;
		/**
		 * @brief operator +
		 * @param other
		 * @return
		 */
        polynomial<T> operator+(const T& other) const;
		/**
		 * @brief operator -
		 * @param other
		 * @return
		 */
        polynomial<T> operator-(const T& other) const;
		/**
		 * @brief operator *
		 * @param other
		 * @return
		 */
        polynomial<T> operator*(const T& other) const;
        /**
         * @brief operator /
         * @param other
         * @return
         */
        polynomial<T> operator/(const T& other) const;

        /******************************/
        /*UNARY OPERATORS             */
        /******************************/
        /**
		 * @brief operator +
		 * @return
		 */
        polynomial<T> operator+() const;
		/**
		 * @brief operator -
		 * @return
		 */
        polynomial<T> operator-() const;

        /******************************/
        /*ASSIGNEMENT (with operators)*/
        /******************************/
		/**
		 * @brief operator =
		 * @param other
		 * @return
		 */
        polynomial<T>& operator=(const polynomial<T> &other);
		/**
		 * @brief operator =
		 * @param other
		 * @return
		 */
        polynomial<T>& operator=(const T &other);

        /**
		 * @brief operator +=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator+=(const polynomial<T> &other);
		/**
		 * @brief operator -=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator-=(const polynomial<T> &other);
		/**
		 * @brief operator *=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator*=(const polynomial<T> &other);
		/**
		 * @brief operator +=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator+=(const T& other);
		/**
		 * @brief operator -=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator-=(const T& other);
		/**
		 * @brief operator *=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator*=(const T& other);
		/**
		 * @brief operator /=
		 * @param other
		 * @return
		 */
        polynomial<T>& operator/=(const T& other);

        /******************************/
        /*COMPARISON                  */
        /******************************/
        /**
		 * @brief operator ==
		 * @param other
		 * @return
		 */
        bool operator==(const polynomial<T> &other) const;
		/**
		 * @brief operator !=
		 * @param other
		 * @return
		 */
        bool operator!=(const polynomial<T> &other) const;

        /******************************/
        /*I/O OPERATOR                */
        /******************************/
        /**
         * @brief operator <<
         * @param os
         * @param poly
         * @return
         */
        friend ostream &operator<<(ostream &os, const polynomial<T> &poly) {
            std::vector<T> coeffs = poly.get_coeffs();
            int nvar = poly.get_nvar();
            int idx=0;
            os <<std::setfill(' ')<<setw(16);
            for(int i=0; i<nvar; i++)
                os << poly.get_basis_name() << "(x"<<i<<")\t";
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
        /******************************/
        /*EVALUATION & COMPOSITION    */
        /******************************/
        /**
         * @brief evaluate_1Dbase
         * @param other
         * @return
         */
        std::vector<polynomial<T> > evaluate_base1D(const polynomial<T> &other) const;
        /**
         * @brief evaluate_1Dbase
         * @param other
         * @return
         */
        std::vector<T> evaluate_base1D(const T &other) const;
        /**
         * @brief composition
         * @param other
         * @return
         */
        polynomial<T> composition(const std::vector<polynomial<T> > &other) const;
        /**
         * @brief evaluate
         * @param x
         * @return
         */
        T evaluate(const std::vector<T> &x) const;
    public:
        //getter and setters
        /**
         * @brief get_name
         * @return
         */
        std::string get_name(){return m_name;} //name of the derived class like "Chebyshev_Polynomial"

        std::string get_basis_name(){return "X";}
        /**
         * @brief get_coeffs
         * @return
         */
        std::vector<T> get_coeffs() const {return m_coeffs;}
        /**
         * @brief set_coeffs
         * @param coeffs
         */
        void set_coeffs(std::vector<T> &coeffs){
                if(m_coeffs.size()!=coeffs.size()){
                    std::cout<<"Coefficients vectors don't have the same lenght"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                m_coeffs=coeffs;
        }
        /**
         * @brief set_coeffs
         * @param idx
         * @param value
         */
        void set_coeffs(const int &idx, const T &value){m_coeffs[idx]=value;}
        /**
         * @brief get_degree
         * @return
         */
        int get_degree() const {return m_degree;}
        /**
         * @brief get_nvar
         * @return
         */
        int get_nvar() const {return m_nvar;}
        /**
         * @brief get_range
         * @return
         */
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
        //Matrix needed for the multivariate polynomial ordering
        void initialize_J();
        void initialize_N();
        //Matrix need for fast multiplication
        void initialize_M(const int &nvar, const int &degree);
        void delete_M();

    protected:
        std::vector<std::vector<int> > get_J() const {return m_J;}
        std::vector<std::vector<int> > get_N() const {return m_N;}
        std::vector<int> get_row(const int &idx, const int &deg) const;
        int get_idx(const std::vector<int> &k) const;

    protected:
        string m_name;
        std::vector<T> m_coeffs;
        int m_degree;
        int m_nvar;

        std::vector<std::vector<int> > m_J, m_N;
		static std::vector<int> m_M;
		static int m_Mnvar, m_Mdegree;
	};
	
	template < class T>
    static polynomial<T> operator-(const T left, const polynomial<T> right){
		return -right+left;
	}

	template < class T>
    static polynomial<T> operator-(const int left, const polynomial<T> right){
		return -right+((T) left);
	}

	template < class T>
    static polynomial<T> operator+(const T left, const polynomial<T> right){
		return right+left;
	}

	template < class T>
    static polynomial<T> operator+(const int left, const polynomial<T> right){
		return right+((T) left);
	}

	template < class T>
    static polynomial<T> operator*(const T left, const polynomial<T> right){
		return right*left;
	}

	template < class T>
    static polynomial<T> operator*(const int left, const polynomial<T> right){
		return right*((T) left);
	}

}}



#endif /* polynomial_H_ */
