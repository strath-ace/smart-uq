#include "main_list.h"

void main_vanderpol()
{
    std::ofstream file;
    file.open ("results.out");

    //algebra params
    int degree = 10;
    int nvar = 2;
    //integration params
    double step = 0.01;
    double tend = 5.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    ranges[0][0] = -2.0; ranges[0][1] = 2.0;
    ranges[1][0] = -2.0; ranges[1][1] = 2.0;

    std::vector<Chebyshev_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Chebyshev_Polynomial<double>(nvar,degree));
        x0[i].set_coeffs(i+1,1);
    }

    std::vector<Chebyshev_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Chebyshev_Polynomial<double>(nvar,degree));
    }

    //translation  [-1,1] ----> [a,b]
    for(int i=0; i<nvar; i++){
        x0[i] = (ranges[i][1]-ranges[i][0])/2.0*x0[i] + (ranges[i][1]+ranges[i][0])/2.0;
    }

    //assign initial status
    for(int i=0; i<nvar; i++){
        res[i] = x0[i];
    }

    //perform integration
    for(int i=0; i<tend/step; i++){
        std::cout<<"iteration "<<i<<std::endl;


        res = euler<double>(f,res,param0,step);

        for(int j=0; j<nvar; j++){
            std::vector<double> coeffs = res[j].get_coeffs();
            for(int k=0; k<res[j].get_coeffs().size(); k++){
                file << coeffs[k] << " ";
            }
            file << "\n";
        }

    }
    file.close();
}

