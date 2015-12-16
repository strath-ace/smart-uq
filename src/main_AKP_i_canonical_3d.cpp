#include "main_list.h"

void main_AKP_i_canonical(){
    std::ofstream file;
    file.open ("intrusive_canonical_case3_3D.txt");
for (int degree = 4; degree <= 4; degree ++){
    //algebra params
    // int degree = 4;
    int nvar = 6;
    int nparam = 3; //*********
    //integration params
    double step = 0.01;
    double sma = 1; //*********
    double tend = 2.0*M_PI/pow(sma,-3.0/2.0);
    double e = 0; //*********

    std::vector<std::vector<double> > ranges_x, ranges_p;
    for(int i=0; i<nvar; i++){
        ranges_x.push_back(std::vector<double>(2));
        ranges_x[i][0] = -1.0; ranges_x[i][1] = 1.0;
    }
    for(int i=0; i<nparam; i++){
        ranges_p.push_back(std::vector<double>(2));
        ranges_p[i][0] = -1.0; ranges_p[i][1] = 1.0;
    }

    std::vector<double> x(nvar), param(nparam), unc_x(nvar), unc_p(nparam);

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = sqrt(1+e);
    x[5] = 0.0;

    if(nparam>0)
        param[0] = 0.01;

    unc_x[0] = 0.01;
    unc_x[1] = 0.01;
    unc_x[2] = 0.01;
    unc_x[3] = 0.005;
    unc_x[4] = 0.005;
    unc_x[5] = 0.005;

    for (int i=0; i<nparam; i++)
        unc_p[i] = param[0]* 10.0/100.0; //10% of uncertainty on the model parameter

    for(int i=0; i<nvar; i++){
        ranges_x[i][0] = x[i]-unc_x[i];
        ranges_x[i][1] = x[i]+unc_x[i];
    }
    for(int i=0; i<nparam; i++){
        ranges_p[i][0] = param[i]-unc_p[i];
        ranges_p[i][1] = param[i]+unc_p[i];
    }

    clock_t begin,end;
    begin=clock();

    std::vector<Canonical_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
        x0[i].set_coeffs(i+1,1);
    }
    for(int i=0; i<nparam; i++){
        param0.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
        param0[i].set_coeffs(i+1+nvar,1);
    }

    std::vector<Canonical_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
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
            // std::cout<<"iteration "<<i<<std::endl;

            res = rk4(f,res,param0,step);

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
            if((k+1)%nvar == 0  && nparam>0){
                std::vector<double> param_coeff = param0[0].get_coeffs();
                for(int kk=0; kk<param_coeff.size(); kk++)
                    file << setprecision(16) << param_coeff[kk] << " ";
                file << "\n";
            }
        }
        file.close();
        exit(EXIT_FAILURE);
    }
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "canonical, time elapsed : " << time_akp << endl << endl;

    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++)
            file << setprecision(16) << coeffs_all[k][kk] << " ";
        file << "\n";
        if((k+1)%nvar == 0 && nparam>0){
            for (int p=0;p<nparam;p++){
                std::vector<double> param_coeff = param0[p].get_coeffs();
                for(int kk=0; kk<param_coeff.size(); kk++)
                    file << setprecision(16) << param_coeff[kk] << " ";
                file << "\n";
            }
        }
    }
file << "\n\n\n\n";
}
    file.close();

}