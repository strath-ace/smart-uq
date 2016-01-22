#include "../include/smartuq.h"

int main()
{
    //algebra params
    int degree = 3;
    int nvar = 2;
    int nparam=0;
    //sampling
    int npoints=combination(nvar+nparam,degree);
    int ncoeffs=combination(nvar+nparam,degree);
    //integration params
    double step = 0.01;
    double tend = 10.0;
    int freq = 10; //every how many iterations we save the results

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar+nparam; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    ranges[0][0] = -2; ranges[0][1] = 2;
    ranges[1][0] = -2; ranges[1][1] = 2;

    //timer
    clock_t begin,end;
    begin=clock();

    //sampling and evaluating basis
    sampling::lhs<double> lhs_gen(nvar+nparam,npoints);
    Eigen::MatrixXd base_matrix (npoints,ncoeffs);

    std::vector<std::vector<double> > x0;
    std::vector<std::vector<double> > param0;

    for(int i=0; i<npoints; i++){
        std::vector<double> xp_aux=lhs_gen();
        std::vector<double> x_next,p_next;
        // translation to real range for x0 and to -1,1 for evaluation of basis
        for (int j=0; j<nvar; j++){
            x_next.push_back(ranges[j][0]+xp_aux[j]*(ranges[j][1]-ranges[j][0]));
            xp_aux[j]*=2;
            xp_aux[j]-=1;
        }
        for (int j=0; j<nparam; j++){
            p_next.push_back(ranges[j+nvar][0]+xp_aux[j+nvar]*(ranges[j+nvar][1]-ranges[j+nvar][0]));
            xp_aux[j+nvar]*=2;
            xp_aux[j+nvar]-=1;
        }
        x0.push_back(x_next);
        param0.push_back(p_next);
        //evaluation of the base
        Chebyshev_Polynomial<double> base_poly(nvar+nparam,degree,1.0);
        base_matrix(i,0)=1.0;
        for (int j=1;j<ncoeffs;j++){
            base_poly.set_coeffs(j-1,0.0);
            base_poly.set_coeffs(j,1.0);
            base_matrix(i,j)=base_poly.evaluate(xp_aux);
        }
    }

    // Eigen::MatrixXd base_inv (npoints,ncoeffs);
    // base_inv=base_matrix.inverse();

    //matrix condition number
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(base_matrix);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    cout<< "CONDITION NUMBER ="<<cond << endl;

    // assign initial status
    std::vector<std::vector<double> > res = x0;

    //perform integration
    std::vector<std::vector<double> > coeffs_all;
    // // //pointwise results
    // for (int v=0; v<nvar; v++){
    //     std::vector<double> coeffs;
    //     for (int j=0;j<npoints;j++){
    //         coeffs.push_back(x0[j][v]);
    //     }
    //     coeffs_all.push_back(coeffs);
    // }

    for(int i=0; i<tend/step; i++){
        // std::cout<<"iteration "<<i<<std::endl;
        for (int j=0;j<npoints;j++){
            res[j] = euler(f,res[j],param0[j],step);
        }
        // //pointwise results
        // for (int v=0; v<nvar; v++){
        //     std::vector<double> coeffs;
        //     for (int j=0;j<npoints;j++){
        //         coeffs.push_back(res[j][v]);
        //     }
        //     coeffs_all.push_back(coeffs);
        // }

        if((i+1)%freq == 0){
            for (int v=0; v<nvar; v++){
                Eigen::VectorXd y(npoints);
                for (int j=0;j<npoints;j++){
                    y(j)=res[j][v];
                }

                Eigen::VectorXd coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
                // Eigen::VectorXd coe = base_inv*y;
                std::vector<double> coeffs;
                for (int c=0;c<ncoeffs;c++){
                    coeffs.push_back(coe(c));
                }
                coeffs_all.push_back(coeffs);
            }
        }
    }
    //timer
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "vanderpol non-intrusive full, time elapsed : " << time_akp << endl << endl;

    // write to file
    std::ofstream file;
    file.open ("vanderpol_ni_euler_n10.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}

