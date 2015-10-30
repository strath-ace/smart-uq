#include "main_list.h"
#include "utils.h"
#include "Sampling/sampling.h"
#include "Eigen/Dense"
#include "f_ni.h"
#include "integrators_ni.h"


void main_AKP_ni(){//coded for nparam=0

    std::ofstream file;
    file.open ("non_intrusive_case1.txt");
    //algebra params
    for (int degree = 4; degree <=12; degree++){

        int nvar = 4;
        int nparam = 0; //*********
        //sampling
        int npoints=combination(nvar,degree);
        int ncoeffs=combination(nvar,degree);
        //integration params
        double step = 0.01;
        double sma = 1.0; //*********
        double tend = 1.01; //2.0*M_PI/pow(sma,-3.0/2.0);
        double e = 0.0; //*********

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
        x[3] = sqrt(1+e);

        if(nparam>0)
            param[0] = 0.01;

        unc_x[0] = 0.01;
        unc_x[1] = 0.01;
        unc_x[2] = 0.005;
        unc_x[3] = 0.005;


        if(nparam>0)
            unc_p[0] = param[0]* 10.0/100.0; //10% of uncertainty on the model parameter

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
        //sampling and evaluating basis
        sampling::lhs<double> lhs_gen(nvar,npoints);
        Eigen::MatrixXd base_matrix (npoints,ncoeffs);

        std::vector<std::vector<double> > x0;
        std::vector<double> param0;
        for(int i=0; i<npoints; i++){
            std::vector<double> x_aux=lhs_gen();
            std::vector<double> x_next;
            // translation to real range for x0 and to -1,1 for evaluation of basis
            for (int j=0; j<nvar; j++){
                x_next.push_back(ranges_x[j][0]+x_aux[j]*(ranges_x[j][1]-ranges_x[j][0]));
                x_aux[j]*=2;
                x_aux[j]-=1;
            }
            x0.push_back(x_next);
            //evaluation of the base
            Chebyshev_Polynomial<double> base_poly(nvar,degree,1.0);
            base_matrix(i,0)=1.0;
            for (int j=1;j<ncoeffs;j++){
                base_poly.set_coeffs(j-1,0.0);
                base_poly.set_coeffs(j,1.0);
                base_matrix(i,j)=base_poly.evaluate(x_aux);
            }
        }


        // for(int i=0; i<nparam; i++){
        //     param0.push_back(Chebyshev_Polynomial<double>(nvar+nparam,degree));
        //     param0[i].set_coeffs(i+1+nvar,1);
        // }
        
        // //translation  [-1,1] ----> [a,b]
        // for(int i=0; i<nvar; i++){
        //     x0[i] = (ranges_x[i][1]-ranges_x[i][0])/2.0*x0[i] + (ranges_x[i][1]+ranges_x[i][0])/2.0;
        // }
        // for(int i=0; i<nparam; i++){
        //     param0[i] = (ranges_p[i][1]-ranges_p[i][0])/2.0*param0[i] + (ranges_p[i][1]+ranges_p[i][0])/2.0;
        // }

        std::vector<std::vector<double> > res = x0;

        // //assign initial status
        // for(int i=0; i<nvar; i++){
        //     res[i] = x0[i];
        // }

        std::vector<std::vector<double> > coeffs_all;
        
        try{
            //perform integration
            for(int i=0; i<tend/step; i++){
                for (int j=0;j<npoints;j++){
                    res[j] = rk4_ni(f_ni,res[j],param0,step);
                }

                if((i+1)%100 == 0){
                    for (int v=0; v<nvar; v++){
                        Eigen::VectorXd y(npoints);
                        for (int j=0;j<npoints;j++){
                            y(j)=res[j][v];
                        }
                        Eigen::VectorXd coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
                        std::vector<double> coeffs;
                        for (int c=0;c<ncoeffs;c++){
                            coeffs.push_back(coe(c));
                        }
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
                // if((k+1)%nvar == 0  && nparam>0){
                //     std::vector<double> param_coeff = param0[0];
                //     for(int kk=0; kk<param_coeff.size(); kk++)
                //         file << setprecision(16) << param_coeff[kk] << " ";
                //     file << "\n";
                // }
            }
            file.close();
            exit(EXIT_FAILURE);
        }

        end=clock();
        double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
        cout << "time elapsed : " << time_akp << endl << endl;
        
        for(int k=0; k<coeffs_all.size(); k++){
            for(int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
            // if((k+1)%nvar == 0 && nparam>0){
            //     std::vector<double> param_coeff = param0[0];
            //     for(int kk=0; kk<param_coeff.size(); kk++)
            //         file << setprecision(16) << param_coeff[kk] << " ";
            //     file << "\n";
            // }
        }

        file << "\n\n\n\n";
    }
    file.close();
}

