#include "main_list.h"
#include <time.h>
#include <ctime>

void main_test_dct(){
    // //INTERMEDIATE OUTPUT
    // std::ofstream file;
    // std::string filename = "intermediate_dcts.out";
    // file.open (filename.data());
    
    // cout << endl << endl << "degree = 4-6   nvar = 1-9" << endl;

    for (int degree=10; degree <=10; degree+=1){
        cout << "[" << endl;
        for (int nvar=4; nvar<=4; nvar++){
            //TIME
            clock_t direct0, directf, dct0, dctf;
            double direct_t, dct_t, t_rel;

            //int nvar = 4;
            //int degree = 10;

            srand(time(NULL));
            int ncoeffs=combination(nvar,degree);

            //initialise polynomials to be multiplied
            std::vector<Chebyshev_Polynomial<double> > x, f;
            for(int i=0; i<2; i++){
                x.push_back(Chebyshev_Polynomial<double>(nvar,degree));
                for (int j=0; j<ncoeffs;j++){
                    x[i].set_coeffs(j,(double) (rand()%1000)/50.0-10); //randomly
                    // x[i].set_coeffs(j,1.0); //with a pattern
                }
            }

            std::vector<double> x0=x[0].get_coeffs();
            std::vector<double> x1=x[1].get_coeffs();

            // // print polynomials with terms detailed
            // cout << x[0] << endl << endl;
            // cout << x[1] << endl;

            // // OUTPUT
            // cout << "____x0____" << endl;
            // for(int k=0; k<x0.size(); k++){
            //     cout <<x0[k] << "   ";
            // }
            // cout << endl << "____x1____" << endl;
            // for(int k=0; k<x1.size(); k++){
            //     cout <<x1[k] << "   ";
            // }
            // cout << endl << endl << "*******************" << endl << "   DIRECT METHOD" << endl << "*******************" << endl;

            // TIME
            direct0=clock();

            Chebyshev_Polynomial<double> x0x1d = x[0]*x[1];
           
            // TIME
            directf=clock();

            std::vector<double> x0x1dcoeffs = x0x1d.get_coeffs();

            // //OUTPUT
            // cout << "_x0*x1_ : ";
            // for(int k=0; k<x0x1dcoeffs.size(); k++){
            //     cout <<x0x1dcoeffs[k] << "   ";
            // }


            // // DCT METHOD
            // cout << endl << endl << "*******************" << endl << "        DCT     " << endl << "*******************" << endl;

            //TIME
            dct0=clock();

            Chebyshev_Polynomial<double> x0x1 = direct_multiplication(x[0],x[1]);

            // TIME
            dctf=clock();

            //get results coefficients
            std::vector<double> x0x1coeffs=x0x1.get_coeffs();

            //OUTPUT
            // cout << endl << "_x0*x1_ : ";
            // for (int i=0;i<ncoeffs;i++){
            //     cout << x0x1[i] << "    ";
            // }
            // cout << endl;

            // //PRECISION WARNING FOR TIME TESTS AND/OR SAVE ERROR
            // bool warned=false;
            std::vector<double> acc_rel;
            double max_err=0;
            for (int i=0;i<ncoeffs;i++){
                if (fabs(x0x1coeffs[i])<=ZERO && fabs(x0x1dcoeffs[i])<=ZERO){acc_rel.push_back(1.0);}
                else{
                    acc_rel.push_back(x0x1coeffs[i]/x0x1dcoeffs[i]);
                }

                // // PRECISION WARNING FOR TIME TESTS
                // if ( warned==false && ( err > pow(10,-8) )) {
                //     cout << endl << endl
                //         << "**********************" << endl
                //         << "**********************" << endl
                //         << "** PRECISION  ISSUE **" << endl
                //         << "**********************" << endl
                //         << "**********************" << endl;
                //     warned=true;
                    
                // }

                // ERROR
                double err = fabs(acc_rel[i]-1.0);
                if ( err > max_err ) max_err = err;
            }

            // //PRECISION TEST
            // cout << endl << endl << "DCT/DIRECT : ";
            // for (int i=0;i<ncoeffs;i++){
            //     cout << acc_rel[i] << "    ";
            // }
            // cout << endl << endl;

            //TIME
            direct_t=(double (directf-direct0))/CLOCKS_PER_SEC;
            dct_t=(double (dctf-dct0))/CLOCKS_PER_SEC;
            t_rel=direct_t/dct_t;
            cout << setprecision(16)<< "[" << direct_t << "," << dct_t << "," << t_rel << "," << max_err <<"],"<<endl;
        }
    cout << "]" << endl;
    }
    // //INTERMEDIATE OUTPUT
    // file.close();
}
