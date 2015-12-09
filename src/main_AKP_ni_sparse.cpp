#include "main_list.h"
// #include <Eigen/SVD>

std::vector <std::vector <double> > sparse_to_full(std::vector< std::vector <double> > coeffs_sparse, std::vector < std::vector <int> > index, int dim, int degree){//should add sanity checks and so...
    std::vector < std::vector <double> > coeffs_full;
    for (int p=0; p<coeffs_sparse.size(); p++){
        std::vector <double> coeffs = coeffs_sparse[p];
        Chebyshev_Polynomial<double> poly(dim,degree);
        for (int i=0;i<coeffs.size();i++){
            int deg = std::accumulate(index[i].begin(),index[i].end(),0);
            int pos0=0;
            if (deg>0) pos0=combination(dim,deg-1);
            int pos= poly.get_idx(index[i]);
            poly.set_coeffs(pos+pos0,coeffs[i]);
        }
        coeffs_full.push_back(poly.get_coeffs());
    }
    return coeffs_full;
}

void main_AKP_ni_sparse(){

    std::ofstream file;
    file.open ("non_intrusive_sparse_case4.txt");
    //algebra params
    for (int level = 2; level <=2; level++){

        int degree=pow(2,level);
        int nvar = 4;
        int nparam = 1; //*********
        //sampling
        int ncoeffs = sparse_grid_cfn_size ( nvar+nparam,level );
        int npoints = ncoeffs;
        //integration params
        double step = 0.01;
        double sma = 2.0; //*********
        double tend = 15.01;//2.0*M_PI/pow(sma,-3.0/2.0);
        double e = 0.5; //*********

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

        unc_x[0] = 0.005;
        unc_x[1] = 0.005;
        unc_x[2] = 0.0001;
        unc_x[3] = 0.0001;


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

        //sampling and evaluating basis
        double *weights, *points;
        weights = new double[npoints];
        points = new double[(nvar+nparam)*npoints];

        sparse_grid_cc (nvar+nparam, level, npoints, weights, points);
        
        vector< vector<int> > idx;
        idx = sparse_grid_index( nvar+nparam, level);

        Eigen::MatrixXd base_matrix (npoints,ncoeffs);

        std::vector<std::vector<double> > x0, param0;

        for(int i=0; i<npoints; i++){

            std::vector<double> xp_aux, x_next, p_next;
            for (int j=0; j<nvar+nparam;j++){
                xp_aux.push_back(points[i*(nvar+nparam)+j]);
            }

            //translation to range
            for (int j=0; j<nvar; j++){
                x_next.push_back((ranges_x[j][0]+ranges_x[j][1])/2.0+xp_aux[j]*(ranges_x[j][1]-ranges_x[j][0])/2.0);
            }
            for (int j=0; j<nparam; j++){
                p_next.push_back((ranges_p[j][0]+ranges_p[j][1])/2.0+xp_aux[j+nvar]*(ranges_p[j][1]-ranges_p[j][0])/2.0);
            }

            //storage
            x0.push_back(x_next);
            param0.push_back(p_next);
            
            //evaluation of the base
            for (int c=0; c<ncoeffs; c++){
                double term=1;
                for (int j=0; j<nvar+nparam; j++){
                    Chebyshev_Polynomial<double> base_poly(1,degree);
                    base_poly.set_coeffs(idx[c][j],1);
                    term*=base_poly.evaluate(xp_aux[j]);
                }
                base_matrix(i,c)=term;
            }
        }

        // //check matrix conditioning
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(base_matrix);
        // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
        // cout<< "CONDITION NUMBER ="<< cond << endl;

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

        std::vector<std::vector<double> > coeffs_all_sparse;
        
        try{
            //perform integration
            for(int i=0; i<tend/step; i++){
                for (int j=0;j<npoints;j++){
                    res[j] = rk4(f,res[j],param0[j],step);
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
                        coeffs_all_sparse.push_back(coeffs);
                    }
                }
            }
        }
        catch(const std::exception&)
        {
            std::vector<std::vector<double> > coeffs_all;
            coeffs_all=sparse_to_full(coeffs_all_sparse,idx,nvar+nparam,degree);
            int ncoeffs_full=coeffs_all[0].size();
            for(int k=0; k<coeffs_all.size(); k++){
                for(int kk=0; kk<coeffs_all[k].size(); kk++)
                    file << setprecision(16) << coeffs_all[k][kk] << " ";
                file << "\n";
                if((k+1)%nvar == 0  && nparam>0){// implement this better... i left it like this for postprocessing needs
                    file << "0.01 0 0 0 0 0.0009999999999999992 ";
                    for (int c=6; c<ncoeffs_full; c++) file << "0 ";
                    file << "\n";  
                    for (int p=1;p<nparam;p++){
                        for (int c=0;c<ncoeffs_full;c++){
                            if (c==nvar+p+1) file << "0.001 ";
                            else file << "0 ";
                        }
                        file << "\n";
                    }
                    // std::vector<double> param_coeff = param0[0];
                    // for(int kk=0; kk<param_coeff.size(); kk++)
                    //     file << setprecision(16) << param_coeff[kk] << " ";
                }
            }
            file.close();
            exit(EXIT_FAILURE);
        }

        end=clock();
        double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
        cout << "non-intrusive sparse, time elapsed : " << time_akp << endl << endl;
        
        std::vector<std::vector<double> > coeffs_all;
        coeffs_all=sparse_to_full(coeffs_all_sparse,idx,nvar+nparam,degree);
        int ncoeffs_full=coeffs_all[0].size();
        for(int k=0; k<coeffs_all.size(); k++){
            for(int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
            if((k+1)%nvar == 0  && nparam>0){// implement this better... i left it like this for postprocessing needs
                file << "0.01 0 0 0 0 0.0009999999999999992 ";
                for (int c=6; c<ncoeffs_full; c++) file << "0 ";
                file << "\n";
                for (int p=1;p<nparam;p++){
                    for (int c=0;c<ncoeffs_full;c++){
                        if (c==nvar+p+1) file << "0.001 ";
                        else file << "0 ";
                    }
                    file << "\n";
                }
                // std::vector<double> param_coeff = param0[0];
                // for(int kk=0; kk<param_coeff.size(); kk++)
                //     file << setprecision(16) << param_coeff[kk] << " ";
            }
        }

        file << "\n\n\n\n";
    }
    file.close();
}

