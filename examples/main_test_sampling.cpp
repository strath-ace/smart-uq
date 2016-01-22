#include "../include/smartuq.h"

int main(){
    int dimension = 10;
    int points = 15;
    cout << "SOBOL SAMPLING - "<< points<<" POINTS IN "<< dimension <<"D" << endl;
    sampling::sobol<double> sobol_gen(dimension);
    for (int i=0;i<points;i++){
        std::vector<double> nextpoint=sobol_gen();
        for (int j=0;j<dimension;j++){
            cout << nextpoint[j] << "   ";
        }
        cout << endl;
    }    

    cout << endl;
    cout << "LHS SAMPLING - "<< points<<" POINTS IN "<< dimension <<"D" << endl;
    sampling::lhs<double> lhs_gen(dimension,points);
    for (int i=0;i<points;i++){
        std::vector<double> nextpoint=lhs_gen();
        for (int j=0;j<dimension;j++){
            cout << nextpoint[j] << "   ";
        }
        cout << endl;
    }
}
