#include "../include/smartuq.h"
#include <time.h>
#include <ctime> 

int main(){

    for (int nvar=1; nvar <=9; nvar+=1){
        cout << "nvar = "<< nvar << "   degree 1-8"<< endl;
        cout << "[" << endl;
        for (int degree=1; degree<=8; degree++){
            //TIME
            clock_t direct0, directf;
            double direct_t;

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
            //TIME
            direct0=clock();

            x0x1 = x[0]*x[1];

            // TIME
            directf=clock();

            direct_t=(double (directf-direct0))/CLOCKS_PER_SEC;
            cout << setprecision(16)<< ncoeffs << " , " <<direct_t << " ,"<<endl;
        }

    cout << "]" << endl;

    }
}
