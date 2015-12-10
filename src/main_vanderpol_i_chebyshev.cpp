#include "main_list.h"

void main_vanderpol_i_chebyshev()
{
    //algebra params
    int degree = 10;
    int nvar = 2;
    int nparam=0;
    //integration params
    double step = 0.01;
    double tend = 5.0;
    int freq = 1; //every how many iterations we save the results

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar+nparam; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    ranges[0][0] = -2.0; ranges[0][1] = 2.0;
    ranges[1][0] = -2.0; ranges[1][1] = 2.0;

    //timer
    clock_t begin,end;
    begin=clock();

    std::vector<Chebyshev_Polynomial<double> > x0, param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        x0[i].set_coeffs(i+1,1);
    }
    for(int i=0; i<nparam; i++){
        param0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        param0[i].set_coeffs(nvar+i+1,1);
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
        param0[i] = (ranges[nvar+i][1]-ranges[nvar+i][0])/2.0*param0[i] + (ranges[nvar+i][1]+ranges[nvar+i][0])/2.0;
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
                // Canonical_Polynomial<double> res_canon(nvar,degree);
                // res_canon.assign_from_chebyshev(coeffs);
                // coeffs_all.push_back(res_canon.get_coeffs());
            }
            // // depends on the postprocessor for nparam!=0
            // for(int j=0; j<nparam; j++){
            //     std::vector<double> coeffs = param0[j].get_coeffs();
            //     coeffs_all.push_back(coeffs);
            // }
        }
    }
    //timer
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "vanderpol chebyshev, time elapsed : " << time_akp << endl << endl;

    // write to file
    std::ofstream file;
    file.open ("results_vanderpol_chebyshev_euler_n10_canonical_t50s.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}

