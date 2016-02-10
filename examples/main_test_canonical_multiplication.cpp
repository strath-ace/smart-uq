#include "../include/smartuq.h"
#include <time.h>
#include <ctime> 

int main(){

    int degmax = 5;
    int nvarmax = 100;

    int degmin = 5;
    int nvarmin = 51;

    std::vector<double> t_D, t_M;
    for (int degree=degmin; degree <=degmax; degree++){
        for (int nvar=nvarmin; nvar<=nvarmax; nvar++){
            
            //TIME
            clock_t begin, end;

            srand(time(NULL));
            int ncoeffs=combination(nvar,degree);

            //initialise polynomials to be multiplied
            std::vector<Canonical_Polynomial<double> > x;
            for(int i=0; i<2; i++){
                x.push_back(Canonical_Polynomial<double>(nvar,degree));
                for (int j=0; j<ncoeffs;j++){
                    x[i].set_coeffs(j,(double) (rand()%1000)/50.0-10); //randomly
                    // x[i].set_coeffs(j,1.0); //with a pattern
                }
            }

            Canonical_Polynomial<double> x0x1_D(nvar,degree);
            Canonical_Polynomial<double> x0x1_M(nvar,degree);

            // //direct
            // begin=clock();
            // x0x1_D = x[0]*x[1];
            // end=clock();

            // t_D.push_back((double (end-begin))/CLOCKS_PER_SEC);

            //M
            Canonical_Polynomial<double>::initialize_M(nvar,degree);
            begin=clock();
            x0x1_M = x[0]*x[1];
            end=clock();
            t_M.push_back((double (end-begin))/CLOCKS_PER_SEC);
            Canonical_Polynomial<double>::delete_M();

            //check equality
            // if (x0x1_M!=x0x1_D) cout << "PRECISION ISSUE @  degree = " << degree <<" ,  nvar = "<< nvar << " !!!" << endl;
        }
    }

    // //printouts
    int count;

    // cout << "t_D = ["<< endl;
    // count = 0;
    // for (int degree=degmin; degree <=degmax; degree++){
    //     cout << "[ ";
    //     for (int nvar=nvarmin; nvar<=nvarmax; nvar++){
    //         cout << setprecision(16) << t_D[count] << (nvar==nvarmax ? " ]" : " ,  ");
    //         count ++;
    //     }
    //     cout << (degree==degmax ? "]" : ",") << endl;
    // }

    cout << "t_M = ["<< endl;
    count = 0;
    for (int degree=degmin; degree <=degmax; degree++){
        cout << "[ ";
        for (int nvar=nvarmin; nvar<=nvarmax; nvar++){
            cout << setprecision(16) << t_M[count] << (nvar==nvarmax ? " ]" : " ,  ");
            count ++;
        }
        cout << (degree==degmax ? "]" : ",") << endl;
    }

    cout << "ratio = ["<< endl;
    count = 0;
    for (int degree=degmin; degree <=degmax; degree++){
        cout << "[ ";
        for (int nvar=nvarmin; nvar<=nvarmax; nvar++){
            cout << setprecision(16) << t_D[count]/t_M[count] << (nvar==nvarmax ? " ]" : " ,  ");
            count ++;
        }
        cout << (degree==degmax ? "]" : ",") << endl;
    }

}
