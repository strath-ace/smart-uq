#include "main_list.h"

void main_AKP(){
    std::ofstream file;
    file.open ("plotAKP.out");

    //algebra params
    int degree = 4;
    int nvar = 4;
    int nparam = 1;

    //integration params
    double step = 0.01;
    double sma = 1;
    double tend = 2.0*M_PI/pow(sma,-3.0/2.0);
    double e = 0.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }
    std::vector<std::vector<double> > ranges_param;
    for(int i=0; i<nparam; i++){
        ranges_param.push_back(std::vector<double>(2));
        ranges_param[i][0] = -1.0; ranges_param[i][1] = 1.0;
    }

    std::vector<double> x(nvar), param(nparam), unc(nvar), unc_param(nparam);

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = sqrt(1+e);

    param[0] = 0.01;

    unc[0] = 0.01;
    unc[1] = 0.01;
    unc[2] = 0.005;
    unc[3] = 0.005;

    unc_param[0] = param[0]* 5.0/10.0;

    ranges[0][0] = x[0]-unc[0];
    ranges[0][1] = x[0]+unc[0];

    ranges[1][0] = x[1]-unc[1];
    ranges[1][1] = x[1]+unc[1];

    ranges[2][0] = x[2]-unc[2];
    ranges[2][1] = x[2]+unc[2];

    ranges[3][0] = x[3]-unc[3];
    ranges[3][1] = x[3]+unc[3];

    ranges_param[0][0] = param[0]-unc_param[0];
    ranges_param[0][1] = param[0]+unc_param[0];

    std::vector<Chebyshev_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        x0[i].set_coeffs(i+1,1);
    }
    for(int i=0; i<nparam; i++){
        param0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        param0[i].set_coeffs(1+nvar+i,1);
    }

    std::vector<Chebyshev_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
    }

    //translation  [-1,1] ----> [a,b]
    for(int i=0; i<nvar; i++){
        x0[i] = (ranges[i][1]-ranges[i][0])/2.0*x0[i] + (ranges[i][1]+ranges[i][0])/2.0;
    }
    for(int i=0; i<nparam; i++){
        param0[i] = (ranges_param[i][1]-ranges_param[i][0])/2.0*param0[i] + (ranges_param[i][1]+ranges_param[i][0])/2.0;
    }

    //assign initial status
    for(int i=0; i<nvar; i++){
        res[i] = x0[i];
    }

    //perform integration
    for(int i=0; i<tend/step; i++){
        std::cout<<"iteration "<<i<<std::endl;

        res = rk4<double>(f,res,param0,step);

        if((i+1)%100 == 0){
            for(int j=0; j<nvar; j++){
                std::vector<double> coeffs = res[j].get_coeffs();
                for(int k=0; k<res[j].get_coeffs().size(); k++){
                    file << setprecision(16) << coeffs[k] << " ";
                }
                file << "\n";
            }
        }
    }
    file.close();
}

