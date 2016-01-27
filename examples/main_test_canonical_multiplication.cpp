#include "../include/smartuq.h"
#include <time.h>
#include <ctime> 

int main(){

    for (int nvar=1; nvar <=50; nvar++){
        cout << "nvar = "<< nvar << "   degree 5"<< endl;
        cout << "[" << endl;
        for (int degree=5; degree<=5; degree++){
            //TIME
            clock_t direct0, directf;
            double direct_t, direct_t1;

            //int nvar = 4;
            //int degree = 10;

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

            Canonical_Polynomial<double> x0x1(nvar,degree);
            Canonical_Polynomial<double> x0x1_1(nvar,degree);

            //TIME
            direct0=clock();

            x0x1 = x[0]*x[1];

            // TIME
            directf=clock();

            direct_t=(double (directf-direct0))/CLOCKS_PER_SEC;
            cout << setprecision(16)<< ncoeffs << " , " <<direct_t << " ,"<<endl;

            Canonical_Polynomial<double>::initialize_M(nvar,degree);
            direct0=clock();
            x0x1_1 = x[0]*x[1];
            directf=clock();
            direct_t1=(double (directf-direct0))/CLOCKS_PER_SEC;
            cout << setprecision(16)<< ncoeffs << " , " <<direct_t1 << " ," << direct_t/direct_t1<<endl;
            cout << (x0x1==x0x1_1 ? "equal" : "f***ing different") << endl;
            Canonical_Polynomial<double>::delete_M();
        }

    cout << "]" << endl;

    }
}
