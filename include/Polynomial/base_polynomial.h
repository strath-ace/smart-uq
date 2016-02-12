/******************************************************************************
 *                           POLYNOMIAL_H                           *
 *            Chebyshev Polynomial Algebra of the SMART-UQ toolbox            *
 ******************************************************************************/

/*
---------------- Copyright (C) 2015 University of Strathclyde----------------
---------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------------
------------------------- Author: Annalisa Riccardi --------------------------
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
#include "utils.h"

#include "exception.h"

using namespace std;

namespace smart{
namespace polynomial{

        /**
         * @brief The base_polynomial class
         */
        template < class T >
	class base_polynomial{
	public:

		static const int MAX_DEGREE = 100;
		
		/**
		 * @brief base_polynomial
		 * @param name
		 * @param vars
		 * @param order
		 */
		base_polynomial(const string &name, const int &vars, const int &order);

		/**
		 * @brief base_polynomial
		 * @param name
		 * @param vars
		 * @param order
		 * @param i
		 */
		base_polynomial(const string &name, const int &vars, const int &order, const int &i);

		/**
		 * @brief base_polynomial
		 * @param name
		 * @param vars
		 * @param order
		 * @param value
		 */
		base_polynomial(const string &name, const int &vars, const int &order, const T &value);

		//deconstructor
		/**
		 * @brief ~base_polynomial
		 */
		virtual ~base_polynomial(){}

		//virtual getters
		/**
		 * @brief get_basis_name
		 * @return
		 */
		virtual std::string get_basis_name() const = 0; //basis identifier like "T" for Chebyshev

		/**
		 * @brief operator <<
		 * @param os
		 * @param poly
		 * @return
		 */
		friend ostream &operator<<(ostream &os, const base_polynomial<T> &poly) {
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
		//getter and setters
		/**
		 * @brief get_coeffs
		 * @return
		 */
		std::vector<T> get_coeffs() const {return m_coeffs;}
		/**
		 * @brief get_name
		 * @return
		 */
		std::string get_name(){return m_name;} //name of the derived class like "Chebyshev_Polynomial"
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
		 * @brief get_idx
		 * @param k
		 * @return
		 */
		int get_idx(const std::vector<int> &k) const;
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

	private:
		//polynomial representation variables
		void initialize_J();
		void initialize_N();

	protected:
		std::vector<std::vector<int> > get_J() const {return m_J;}
		std::vector<std::vector<int> > get_N() const {return m_N;}
		std::vector<int> get_row(const int &idx, const int &deg) const;

	protected:
		string m_name;
		std::vector<T> m_coeffs;
		int m_degree;
		int m_nvar;
		std::vector<std::vector<int> > m_J, m_N;

	};
}
}



#endif /* BASE_POLYNOMIAL_H_ */
