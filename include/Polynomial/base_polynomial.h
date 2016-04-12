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
#include "../exception.h"
#include "../constants.h"
#include "../wrapper.h"
#include "../Math/smart_math.h"

using namespace std;

namespace smartuq{
namespace polynomial{

    /**
     *@brief The base_polynomial class is a template abstract class. Any new polynomial base added to the toolbox needs to inherit from it and implement its virtual methods.
     *
     * The base_polynomial class is a template abstract class. Any polynomial added to the toolbox needs to inherit from it and implement its virtual methods,
     * namely a method to evaluate the corresponding polynomial basis in scalar and vectorial values, a method for mapping polynomial variables in differnet intervals and
     * a method for changing to and from monomial basis (to be used for fast multiplication in polynomial base)
     */
	template < class T >
    class base_polynomial{

	public:
        /**************/
        /*CONSTRUCTORS*/
        /**************/
        /**
         * @brief base_polynomial constructor
         *
         * The constructor initialize a null polynomial for a given number of variables and degree.
         * It is possible through the constructor to define the range of each variable. If range is not specified [-1,1]
         * is assumed
         * @param vars number of variables
         * @param order polynomial maximum degree
         * @param a vector containing the lower bound for each variable (default -1)
         * @param b vector containing the upper bound for each variable (default 1)
         */
        base_polynomial(const int &vars, const int &order, const std::vector<T> &a=std::vector<T>(), const std::vector<T> &b=std::vector<T>());
        /**
         * @brief base_polynomial constructor
         *
         * The constructor initialize the first order polynomial degree in the i-th variable.
         * The coefficents for degree up to the specified order are allocated and set to zero.
         * @param vars number of variables
         * @param order polynomial maximum degree
         * @param i index of first degree variable
         * @param a variable lower bound (default value -1)
         * @param b variable upper bound (default vale 1)
         */
        base_polynomial(const int &vars, const int &order, const int &i, const T &a=-1.0, const T &b=1.0);
        /**
         * @brief base_polynomial constructor
         *
         * The constructor initialize the zero order polynomial degree with given value.
         * The coefficients for degree up to the specified order are allocated and set to zero
         * @param vars number of variables
         * @param order polynomialmaximum degree
         * @param value polynomial constant term
         */
        base_polynomial(const int &vars, const int &order, const T &value);

        /**
         * @brief ~base_polynomial deconstructor
         */
        virtual ~base_polynomial();

        /******************************/
        /*I/O OPERATOR                */
        /******************************/
        /**
         * @brief operator << I/O operator
         *
         * Writing polynomial to std output
         * @param os output stream
         * @param poly polynomial to be printed
         * @return output stream
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
         * @brief evaluate_basis multivariate basis evaluation
         *
         * Evaluate the multivariate polynomial base in the given point
         * @param x point for evaluation
         * @return vector containing the basis evaluated (ordered according to the lexicographic order)
         */
        virtual std::vector<T> evaluate_basis(const std::vector<T> &x) const = 0;
        /**
         * @brief evaluate multivariate polynomial evaluation
         *
         * Evaluate the polynomial in the given point
         * @param x point for evaluation
         * @return polynomial value
         */
        virtual T evaluate(const std::vector<T> &x) const = 0;
        /**
         * @brief evaluate univariate polynomial evaluation
         *
         * Evaluate the univariate polynomial in a given point
         * @param x value for evaluation
         * @return polynomial value
         */
        virtual T evaluate(const T &x) const = 0;

        /******************************/
        /* CHANGE TO MONOMIAL BASIS   */
        /******************************/
        /**
         * @brief to_monomial_basis transformation to monomial base
         *
         * The method transforms the current polynomial into the monomial base.
         * Virtual method.
         */
        virtual void to_monomial_basis() = 0;
        /**
         * @brief from_monomial_basis transformation from monomial base
         *
         * The method transforms the current polynomial from the monomial base to its native one.
         * Virtual method
         */
        virtual void from_monomial_basis() = 0;

        /******************************/
        /*INTERPOLATION               */
        /******************************/
        /**
         * @brief interpolation find the polynomial coefficients by mean of interpolation
         *
         * Given a set of point and a corresponding responses, the coefficients of the current polynomial are set by solving
         * the linear system Hc=y, where H is the matrix containing in each row the basis evaluated in one of the points and y the corresponding
         * responses. If the matrix H is square the system is solved by matrix inversion. If the matrix has more rows than columns the system is
         * solved by a Least Square approach. In both cases the inverse of the matrix H is returned by the method.
         *
         * @param[in] x set of points
         * @param[in] y responses
         * @param[out] Hinv matrix inverse
         */

        void interpolation(const std::vector<std::vector<T> > &x, const std::vector<T>  &y, std::vector<std::vector<T> > &Hinv);

        /**
         * @brief interpolation find the polynomial coefficients by mean of interpolation
         *
         * Same as the above method with the differnece that y is a matrix of N responses. Hence one matrix inversion and N matrix vector
         * multiplication need to be performed. The coeffients of the current polynomial are unchanged.
         * The coeffients of the N polynomials, solution of the N linear systems are returned by the method.
         *
         * TODO: tobe changed to static function
         *
         * @param[in] x set of points
         * @param[in] y set of responses
         * @param[out] Hinv inverse of matrix H
         * @param[out] res_coeffs coeffients of the polynomials correpsonding to the N responses
         */

        void interpolation(const std::vector<std::vector<T> > &x, const std::vector<std::vector<T> >  &y, std::vector<std::vector<T> > &Hinv, std::vector<std::vector<T> > &res_coeffs) const;

        /**
         * @brief solve given the inverse of the matrix H and the matrix of responses return the corresponding vector of coefficients
         *
         * Same as above but doesn't need to compute the inverse of the matrix H. It is given as an input
         * @param[in] Hinv inverse of matrix H
         * @param[in] y set of responses
         * @param[out] res_coeffs coeffients of the polynomials correpsonding to the N responses
         */
        void solve(const std::vector<std::vector<T> > &Hinv, const std::vector<std::vector<T> >  &y, std::vector<std::vector<T> > &res_coeffs) const;

        /******************************/
        /*MAPPING                     */
        /******************************/
        /**
         * @brief map function for mapping the variables of the polynomial from [-1,1] to [a,b]
         * @param a vector defining the lower bound variables
         * @param b vector defining the upper bound variables
         */
        virtual void map(const std::vector<T> &a, const std::vector<T> &b) = 0;

    public:
        //getter and setters
        /**
         * @brief get_name return polynomial name
         * @return na,e
         */
        std::string get_name(){return m_name;} //name of the derived class
        /**
         * @brief get_basis_name return a string that identifies the polynomial base
         *
         * It is used when polynomial is printed to std output
         * @return base identifier
         */
        virtual std::string get_basis_name() const = 0;
        /**
         * @brief get_coeffs return polynomial coefficients
         * @return vector of coefficients
         */
        std::vector<T> get_coeffs() const {return m_coeffs;}
        /**
         * @brief set_coeffs set polynomial coefficients
         * @param coeffs vector of polynomial coefficients
         */
        void set_coeffs(std::vector<T> &coeffs){
                if(m_coeffs.size()!=coeffs.size()){
                    std::cout<<"Coefficients vectors don't have the same lenght"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                m_coeffs=coeffs;
        }
        /**
         * @brief set_coeffs set a coefficient to a specific value
         * @param idx index of the coefficient
         * @param value constant value to be set
         */
        void set_coeffs(const int &idx, const T &value){
            m_coeffs[idx]=value;
        }
        /**
         * @brief get_degree return polynomial degree
         * @return polynomial degree
         */
        int get_degree() const {return m_degree;}
        /**
         * @brief get_nvar return polynomial number of variables
         * @return number of variables
         */
        int get_nvar() const {return m_nvar;}
        /**
         * @brief get_range estimate range of the polynomial
         * @return 2D vector containing the polynomial bounds
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
         * @brief is_monomial_base check if the polynomial has been transformed into monomial base or not
         * @return boolean value
         */
        bool is_monomial_base() const {return m_monomial_base;}

        //Matrix need for fast monomial multiplication
        /**
         * @brief initialize_M initialize matrix for fast polynomial multiplication in monomial base
         *
         * TODO: make it static
         * @param nvar number of polynomial variables
         * @param degree polynomial degree
         */
        void initialize_M(const int &nvar, const int &degree);
        /**
         * @brief delete_M free memory allocated by M
         */
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
         * @brief monomial_multiplication routine for fast polynomial multiplication in monomial base
         * @param[in] x1 first factor of the multiplication
         * @param[in] x2 second factor of the multiplication
         * @param[out] res_poly result of the polynomial multilpication between x1 and x2
         */
        static void monomial_multiplication(const base_polynomial<T> &x1, const base_polynomial<T> &x2, base_polynomial<T> &res_poly);

        /******************************/
        /* EVALUATION IN MONOMIAL     */
        /******************************/
        //multivariate for reals
        /**
         * @brief evaluate_basis_monomial evaluate the multivariate monomial base in the corresponding values
         * @param x point for evaluation
         * @return vector of the evaluation of the monomial base
         */
        std::vector<T> evaluate_basis_monomial(const std::vector<T> &x) const;

        //univariate for reals
        /**
         * @brief horner fast implementation of the evaluation of the i-th term of the univariate polinomial
         * @param x point for evaluation
         * @param i index
         * @return value of the evaluation of the i-th term
         */
        T horner(T x, int i) const;

        /**
         * @brief get_J return indexing matrix for reconstruct polynomial orderning
         * @return matrix J
         */
        std::vector<std::vector<int> > get_J() const {return m_J;}
        /**
         * @brief get_N return indexing matrix for reconstruct polynomial ordering
         * @return matrix N
         */
        std::vector<std::vector<int> > get_N() const {return m_N;}
        /**
         * @brief get_row  internal routines for accessing polynomial terms within the lexicographic oder
         * @param idx
         * @param deg
         * @return
         */
        std::vector<int> get_row(const int &idx, const int &deg) const;
        /**
         * @brief get_idx internal routines for accessing polynomial terms within the lexicographic oder
         * @param k
         * @return
         */
        int get_idx(const std::vector<int> &k) const;

    protected:
        /**
         * @brief m_name polynomial name
         */
        string m_name;
        /**
         * @brief m_coeffsvector of coefficients
         */
        std::vector<T> m_coeffs;
        /**
         * @brief m_degree polynomial degree
         */
        int m_degree;
        /**
         * @brief m_nvar number of variables
         */
        int m_nvar;
        /**
         * @brief m_monomial_base flag for transformation into monomial base
         */
        mutable bool m_monomial_base;

        /**
         * @brief m_J matrix for lexicographic ordering
         */
        std::vector<std::vector<int> > m_J;
        /**
         * @brief m_N matrix for lexicographic ordering
         */
        std::vector<std::vector<int> > m_N;

        /**
         * @brief m_a lower variables bounds
         */
        std::vector<T> m_a;
        /**
         * @brief m_b upper variables bounds
         */
        std::vector<T> m_b;

        /**
         * @brief m_M indexing matrix for fast multiplication
         */
        static std::vector<int> m_M;
        /**
         * @brief m_Mnvar matrix M corresponding number of variables
         */
        static int m_Mnvar;
        /**
         * @brief m_Mdegree matrix M corresponding degree
         */
        static int m_Mdegree;
    };

}}



#endif /* base_polynomial_H_ */
