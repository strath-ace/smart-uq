#include "../include/smartuq.h"

int main()
{
    //algebra params
    int degree = 4;
    int nvar = 2;
    int nparam = 4;
    //integration params
    double step = 0.01;
    double tend = 40.0;
    int freq = 10; //every how many iterations we save the results

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
    x[1] = 0.5;
    param[0] = 1.0;
    param[1] = 1.0;
    param[2] = 1.0;
    param[3] = 1.0;
    
    for (int i=0; i<nvar; i++)
        unc_x[i] = x[i]* 0.05; //uncertainty on the model states

    for (int i=0; i<nparam; i++)
        unc_p[i] = param[i]* 0.05; //uncertainty on the model parameter

    for(int i=0; i<nvar; i++){
        ranges_x[i][0] = x[i]-unc_x[i];
        ranges_x[i][1] = x[i]+unc_x[i];
    }

    for(int i=0; i<nparam; i++){
        ranges_p[i][0] = param[i]-unc_p[i];
        ranges_p[i][1] = param[i]+unc_p[i];
    }

    //timer
    clock_t begin,end;
    begin=clock();

    std::vector<Canonical_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
        x0[i].set_coeffs(i+1,1);
    }
    for(int i=0; i<nparam; i++){
        param0.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
        param0[i].set_coeffs(nvar+i+1,1);
    }

    std::vector<Canonical_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Canonical_Polynomial<double>(nvar+nparam,degree));
    }

    //translation  [-1,1] ----> [a,b]
    for(int i=0; i<nvar; i++){
        x0[i] = (ranges_x[i][1]-ranges_x[i][0])/2.0*x0[i] + (ranges_x[i][1]+ranges_x[i][0])/2.0;
        // x0[i]*=2;
    }
    for(int i=0; i<nparam; i++){
        param0[i] = (ranges_p[i][1]-ranges_p[i][0])/2.0*param0[i] + (ranges_p[i][1]+ranges_p[i][0])/2.0;
        // param0[nvar+i]*=2;
    }

    //assign initial status
    for(int i=0; i<nvar; i++){
        res[i] = x0[i];
    }

    //perform integration
    std::vector< std::vector<double> > coeffs_all; 
    for(int i=0; i<tend/step; i++){
        // std::cout<<"iteration "<<i<<std::endl;

        res = euler(f,res,param0,step);

        if((i+1)%freq == 0){
            for(int j=0; j<nvar; j++){
                std::vector<double> coeffs = res[j].get_coeffs();
                coeffs_all.push_back(coeffs);
            }
            // //depends on the postprocessor for nparam!=0
            // for(int j=0; j<nparam; j++){
            //     std::vector<double> coeffs = param0[j].get_coeffs();
            //     coeffs_all.push_back(coeffs);
            // }
        }
    }
    //timer
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "lotka-volterra intr. canonical, time elapsed : " << time_akp << endl << endl;

    // write to file
    std::ofstream file;
    file.open ("lotka_volterra_canonical.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}

