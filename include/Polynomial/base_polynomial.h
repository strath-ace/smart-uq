/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef BASE_POLYNOMIAL_H_
#define BASE_POLYNOMIAL_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <numeric>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "exception.h"
#include "constants.h"
#include "wrapper.h"
#include "../Math/smart_math.h"

using namespace std;

namespace smart{
namespace polynomial{
	template < class T >
    class base_polynomial{

	public:
        /**************/
        /*CONSTRUCTORS*/
        /**************/
        /**
         * @brief base_polynomial
         * @param vars
         * @param order
         */
        base_polynomial(const int &vars, const int &order, const std::vector<T> &a=std::vector<T>(), const std::vector<T> &b=std::vector<T>());
        /**
         * @brief base_polynomial
         * @param vars
         * @param order
         * @param i
         */
        base_polynomial(const int &vars, const int &order, const int &i);
        /**
         * @brief base_polynomial
         * @param vars
         * @param order
         * @param value
         */
        base_polynomial(const int &vars, const int &order, const T &value);
        /**
         * @brief base_polynomial
         * @param vars
         * @param order
         * @param i
         * @param a
         * @param b
         */
        base_polynomial(const int &vars, const int &order, const int &i, const T &a, const T &b);

        virtual ~base_polynomial();

        /******************************/
        /*I/O OPERATOR                */
        /******************************/
        /**
         * @brief operator <<
         * @param os
         * @param poly
         * @return
         */
        friend ostream &operator<<(ostream &os, const base_polynomial<T> &poly) {

            if(poly.is_monomial_base())
                poly.from_monomial_basis();

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
        /*EVALUATION                  */
        /******************************/
        /**
         * @brief evaluate
         * @param x
         * @return
         */
        virtual std::vector<T> evaluate_basis(const std::vector<T> &x) const = 0;
        /**
         * @brief evaluate
         * @param x
         * @return
         */
        virtual T evaluate(const std::vector<T> &x) const = 0;
        /**
         * @brief evaluate
         * @param x
         * @return
         */
        virtual T evaluate(const T &x) const = 0;

        /******************************/
        /* CHANGE TO MONOMIAL BASIS   */
        /******************************/
        /**
         * @brief to_monomial_basis
         * @return
         */
        virtual void to_monomial_basis() = 0;
        /**
         * @brief from_monomial_basis
         * @param coeffs
         */
        virtual void from_monomial_basis() = 0;

        /******************************/
        /*INTERPOLATION               */
        /******************************/
        /**
         * @brief interpolation TODO: Least Square is solved for every Y vector. No efficient implementation of saving internal Eigne matrices has been implemented
         * @param x
         * @param y
         */
        void interpolation(const std::vector<std::vector<T> > &x, const std::vector<T>  &y);

        /**
         * @brief interpolation
         * @param x
         * @param y
         */
        void interpolation(const std::vector<std::vector<T> > &x, const std::vector<std::vector<T> >  &y, std::vector<std::vector<T> > &res_coeffs) const;

        /******************************/
        /*MAPPING                     */
        /******************************/
        virtual void map(const int &idx, const std::vector<T> &a, const std::vector<T> &b) = 0;

    public:
        //getter and setters
        /**
         * @brief get_name
         * @return
         */
        std::string get_name(){return m_name;} //name of the derived class
        /**
         * @brief get_basis_name
         * @return
         */
        virtual std::string get_basis_name() const = 0;
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
        void set_coeffs(const int &idx, const T &value){
            m_coeffs[idx]=value;
        }
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
        /**
         * @brief is_monomial_base
         * @return
         */
        bool is_monomial_base() const {return m_monomial_base;}

        //Matrix need for fast monomial multiplication
        void initialize_M(const int &nvar, const int &degree);
        void delete_M();

    private:
        //Matrix needed for the multivariate polynomial ordering
        void initialize_J();
        void initialize_N();

    protected:
        /******************************************/
        /*ARITHMETIC OPERATIONS subroutines (*)   */
        /******************************************/
        /**
         * @brief monomial_multiplication
         * @param poly_coeffs
         * @param other_coeffs
         * @param res_coeffs
         */
        static void monomial_multiplication(const base_polynomial<T> &x1, const base_polynomial<T> &x2, base_polynomial<T> &res_poly);

        /******************************/
        /* EVALUATION IN MONOMIAL     */
        /******************************/
        // static void evaluate_base1D_monomial(const int &index, const base_polynomial<T> &other, base_polynomial<T>  &out);

        //multivariate for reals
        std::vector<T> evaluate_basis_monomial(const std::vector<T> &x) const;

        //univariate for reals
        T horner(T x, int i) const;

        std::vector<std::vector<int> > get_J() const {return m_J;}
        std::vector<std::vector<int> > get_N() const {return m_N;}
        std::vector<int> get_row(const int &idx, const int &deg) const;
        int get_idx(const std::vector<int> &k) const;

    protected:
        string m_name;
        std::vector<T> m_coeffs;
        int m_degree;
        int m_nvar;
        mutable bool m_monomial_base;

        std::vector<std::vector<int> > m_J, m_N;
        std::vector<T> m_a, m_b;

        static std::vector<int> m_M;
        static int m_Mnvar, m_Mdegree;
	};

}}



#endif /* base_polynomial_H_ */
