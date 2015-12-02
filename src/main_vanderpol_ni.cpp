#include "main_list.h"

void main_vanderpol_ni()
{
    //algebra params
    int degree = 10;
    int nvar = 2;
    int nparam=0;
    //sampling
    int npoints=combination(nvar+nparam,degree);
    int ncoeffs=combination(nvar+nparam,degree);
    //integration params
    double step = 0.01;
    double tend = 5.0;

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
    
    // assign initial status
    std::vector<std::vector<double> > res = x0;

    //perform integration
    std::vector<std::vector<double> > coeffs_all;
    for(int i=0; i<tend/step; i++){
        // std::cout<<"iteration "<<i<<std::endl;
        for (int j=0;j<npoints;j++){
            res[j] = euler(f,res[j],param0[j],step);
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
    //timer
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "vanderpol non-intrusive full, time elapsed : " << time_akp << endl << endl;

    // write to file
    std::ofstream file;
    file.open ("results_vanderpol_ni.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}

