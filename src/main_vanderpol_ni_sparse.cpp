#include "main_list.h"

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

void main_vanderpol_ni_sparse()
{
    //algebra params
    int level = 3;
    int degree = pow(2,level);
    int nvar = 2;
    int nparam=0;
    //sampling
    int ncoeffs = sparse_grid_cfn_size ( nvar+nparam,level );
    int npoints = ncoeffs;
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
            x_next.push_back((ranges[j][0]+ranges[j][1])/2.0+xp_aux[j]*(ranges[j][1]-ranges[j][0])/2.0);
        }
        for (int j=0; j<nparam; j++){
            p_next.push_back((ranges[j+nvar][0]+ranges[j+nvar][1])/2.0+xp_aux[j+nvar]*(ranges[j+nvar][1]-ranges[j+nvar][0])/2.0);
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

    // assign initial status
    std::vector<std::vector<double> > res = x0;

    //perform integration
    std::vector<std::vector<double> > coeffs_all_sparse;
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
                coeffs_all_sparse.push_back(coeffs);
            }
        }
    }
    //timer
    end=clock();
    double time_akp = (double (end-begin))/CLOCKS_PER_SEC;
    cout << "vanderpol non-intrusive sparse, time elapsed : " << time_akp << endl << endl;

    // write to file
    std::vector<std::vector<double> > coeffs_all;
    coeffs_all=sparse_to_full(coeffs_all_sparse,idx,nvar+nparam,degree);

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

