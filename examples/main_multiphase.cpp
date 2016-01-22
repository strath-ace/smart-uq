#include "../include/smartuq.h"

int main(){
    std::ofstream file;
    file.open ("intrusive_case2_multiphase.out");

    //algebra params
    int degree = 4;
    int nvar = 4;
    int nparam = 1;
    int nphases = 3; // number of manoeauvres
    double ttransfer = 2; //time in between manoeuvres

    //integration params
    double step = 0.01;
    double sma = 1;
    double tend = 2.0*M_PI/pow(sma,-3.0/2.0);
    double e = 0.0;

    std::vector<std::vector<double> > ranges_x, ranges_p;
    for(int i=0; i<nvar; i++){
        ranges_x.push_back(std::vector<double>(2));
    }
    for(int i=0; i<nparam; i++){
        ranges_p.push_back(std::vector<double>(2));
    }

    std::vector<double> x(nvar), param(nparam);
    std::vector<double> unc_x(nvar), unc_p(nparam);

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = sqrt(1+e);

    param[0] = 0.5; // deflection of 0.5 rad

    unc_x[0] = 0.01;
    unc_x[1] = 0.01;
    unc_x[2] = 0.005;
    unc_x[3] = 0.005;

    unc_p[0] = param[0]*1.0/10.0;

    for(int i=0; i<nvar; i++){
        ranges_x[i][0] = x[i]-unc_x[i];
        ranges_x[i][1] = x[i]+unc_x[i];
    }
    for(int i=0; i<nparam; i++){
        ranges_p[i][0] = x[i]-unc_p[i];
        ranges_p[i][1] = x[i]+unc_p[i];
    }

    std::vector<Chebyshev_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        x0[i].set_coeffs(i+1,1);
    }

    param0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
    param0[0].set_coeffs(nvar+1,1);

    std::vector<Chebyshev_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
    }

    //translation  [-1,1] ----> [a,b]
    for(int i=0; i<nvar; i++){
        x0[i] = (ranges_x[i][1]-ranges_x[i][0])/2.0*x0[i] + (ranges_x[i][1]+ranges_x[i][0])/2.0;
    }
    for(int i=0; i<nparam; i++){
        param0[i] = (ranges_p[i][1]-ranges_p[i][0])/2.0*param0[i] + (ranges_p[i][1]+ranges_p[i][0])/2.0;
    }

    //assign initial status
    for(int i=0; i<nvar; i++){
        res[i] = x0[i];
    }

    std::vector<std::vector<double> > coeffs_all;
    try{
        //perform integration
        for(int i=0; i<tend/step; i++){
            std::cout<<"iteration "<<i<<std::endl;

            res = rk4<double>(f,res,param0,step);

            //save values to be printed
            if((i+1)%100 == 0){
                for(int j=0; j<nvar; j++){
                    std::vector<double> coeffs = res[j].get_coeffs();
                    coeffs_all.push_back(coeffs);
                }
            }
        }
    }
    catch(const std::exception&)
    {
        for(int k=0; k<coeffs_all.size(); k++){
            for(int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
            if((k+1)%nvar == 0){
                std::vector<double> param_coeff = param0[0].get_coeffs();
                for(int kk=0; kk<param_coeff.size(); kk++)
                    file << setprecision(16) << param_coeff[kk] << " ";
                file << "\n";
            }
        }
        file.close();
        exit(EXIT_FAILURE);
    }

    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++)
            file << setprecision(16) << coeffs_all[k][kk] << " ";
        file << "\n";
        if((k+1)%nvar == 0){
            std::vector<double> param_coeff = param0[0].get_coeffs();
            for(int kk=0; kk<param_coeff.size(); kk++)
                file << setprecision(16) << param_coeff[kk] << " ";
            file << "\n";
        }
    }
    file.close();

}


