
#include "main_list.h"
#include <Eigen/SVD>


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

Eigen::MatrixXd derivxv(Eigen::MatrixXd xv, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
    Eigen::MatrixXd dxv = A*(B*xv);
    return dxv;
}

Eigen::MatrixXd propagate_euler (Eigen::MatrixXd xv, double step, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
    Eigen::MatrixXd xv_next = xv + step*derivxv(xv,A,B);
    return xv_next;
}

void main_spring_mass_ni_sparse()
{
    // Coded only for uncertainty in states (x and v). It can be in all or just in some of them.

// INPUT

// Algebra
    int level = 2;

// System
    // n_DOF of the system
    int n = 10;
    //int repr[] = {1,6}; //masses to obtain polynomial representations of, NOT IMPLEMENTED

    // Nominal initial conditions
    std::vector<double> x0(n,0.0);
    x0[0]=-.25; x0[n-1]=.25;
    std::vector<double> v0(n,0.0);
    
    // Uncertainty in initial conditions
    std::vector<double> unc_x0(n,0.10);
    std::vector<double> unc_v0(n,0.0);
    // Nominal parameters
    std::vector<double> mass(n,1.0); // mass
    std::vector<double> ks(n+1,1.0); // rigidity
    std::vector<double> cs(n+1,0.0); // damping
    // std::vector<double> ls(n+1,1.0); // natural longitude of the spring
    // std::vector<double> fs(n,0.0); // amplitudes of the exciting forces NOT IMPLEMENTED

    // // Uncertainty in parameters NOT IMPLEMENTED
    // std::vector<double> unc_mass(n,0.10);
    // std::vector<double> unc_ks(n+1,0.0);
    // std::vector<double> unc_cs(n+1,0.0);
    // std::vector<double> unc_fs(n,0.0);    

// Simulation
    double step = 0.005;
    double tend = 50.0;
    int freq = 2500; //every how many iterations we save the results



// INITIALISATIONS

    int degree=pow(2,level);
    std::vector < std::vector <int> > idx(2);
    std::vector <int> dim(2,0);
    int dim_unc = 0;
    for (int i=0;i<n;i++){
        if (fabs(unc_x0[i])>ZERO) {
            idx[0].push_back(i);
            dim[0]++;
            dim_unc++;
        }
        if (fabs(unc_v0[i])>ZERO) {
            idx[1].push_back(i);
            dim[1]++;
            dim_unc++;
        }

        // if (fabs(unc_mass[i])>ZERO) {
        //     idx[2].push_back(i);
        //     dim[2]++;
        //     dim_unc++;
        // }
        // if (fabs(unc_fs[i])>ZERO) {
        //     idx[5].push_back(i);
        //     dim[5]++;
        //     dim_unc++;
        //}
    }

    // for (int i=0;i<n+1;i++){
    //     if (fabs(unc_ks[i])>ZERO) {
    //         idx[3].push_back(i);
    //         dim[3]++;
    //         dim_unc++;
    //     }
    //     if (fabs(unc_cs[i])>ZERO) {
    //         idx[4].push_back(i);
    //         dim[4]++;
    //         dim_unc++;
    //     }
    // }

    if(dim_unc==0){
        std::cout<<"This is a deterministic problem, set at least one uncertain variable"<<std::endl;
        exit(EXIT_FAILURE);
    }

    // Number of points and coefficients
    int ncoeffs = sparse_grid_cfn_size ( dim_unc,level );
    int npoints = ncoeffs;
    
    //Creation of vector of ranges for sampling
    std::vector<std::vector<double> > ranges;
    for (int i=0; i<dim[0];i++){
        std::vector<double> r(2);
        r[0] = x0[idx[0][i]]-unc_x0[idx[0][i]];
        r[1] = x0[idx[0][i]]+unc_x0[idx[0][i]];
        ranges.push_back(r);
    }

    for (int i=0; i<dim[1];i++){
        std::vector<double> r(2);
        r[0] = v0[idx[1][i]]-unc_v0[idx[1][i]];
        r[1] = v0[idx[1][i]]+unc_v0[idx[1][i]];
        ranges.push_back(r);
    }

    //timer
    clock_t begin,end;
    begin=clock();

    //sampling and evaluating basis, there will be 2n polynomial representations of dimension dim_unc
    double *weights, *points;
    weights = new double[npoints];
    points = new double[dim_unc*npoints];

    sparse_grid_cc (dim_unc, level, npoints, weights, points);
    
    vector< vector<int> > idx_sparse;
    idx_sparse = sparse_grid_index( dim_unc , level);

    Eigen::MatrixXd base_matrix (npoints,ncoeffs);

    std::vector<std::vector<double> > par;

    for(int i=0; i<npoints; i++){
        std::vector<double> par_aux, par_next;

        // fill par_aux
        for (int j=0; j<dim_unc;j++){
            par_aux.push_back(points[i*dim_unc+j]);
        }

        // translation to real range
        for (int j=0; j<dim_unc; j++){
            par_next.push_back((ranges[j][0]+ranges[j][1])/2.0+par_aux[j]*(ranges[j][1]-ranges[j][0])/2.0);
        }

        // storage
        par.push_back(par_next);
        
        //evaluation of the base
        for (int c=0; c<ncoeffs; c++){
            double term=1.0;
            for (int j=0; j<dim_unc; j++){
                Chebyshev_Polynomial<double> base_poly(1,degree);
                base_poly.set_coeffs(idx_sparse[c][j],1.0);
                term*=base_poly.evaluate(par_aux[j]);
            }
            base_matrix(i,c)=term;
        }
    }

    // interpolate by inversion
    Eigen::MatrixXd base_inv (npoints,ncoeffs);
    base_inv = base_matrix.inverse();

    // //matrix condition number
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(base_matrix);
    // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    // cout<< "CONDITION NUMBER ="<<cond << endl;

    // assign initial status
    Eigen::MatrixXd xv(2*n,npoints);
    std::vector<int> count(2,0);
    for (int i=0; i<n; i++){

        // x
        if (!idx[0].empty() && idx[0][count[0]]==i){
            for (int j=0; j<npoints; j++){
                xv(i,j)=par[j][count[0]];
                ////DEBUG
                // std::cout<<par[j][count[0]]<<"  "<<xv(i,j)<<endl;
            }
            count[0]++;
        }
        else{
            for (int j=0; j<npoints; j++){
                xv(i,j)=x0[i];
            }
        }

        // v
        if (!idx[1].empty() && idx[1][count[1]]==i){
            for (int j=0; j<npoints; j++){
                xv(n+i,j)=par[j][dim[0]+count[1]];
                ////DEBUG
                // std::cout<<par[j][dim[0]+count[1]]<<"   "<<xv(n+i,j)<<endl;
            }
            count[1]++;
        }
        else{
            for (int j=0; j<npoints; j++){
                xv(n+i,j)=v0[i];
            }
        }
    }

    // build system matrices
    Eigen::SparseMatrix<double> A(2*n,2*n);
    Eigen::SparseMatrix<double> B(2*n,2*n);
    A.reserve(Eigen::VectorXi::Constant(2*n,1));
    B.reserve(Eigen::VectorXi::Constant(2*n,4));

    for (int i=0;i<n;i++){
        int j=n+i;

        A.insert(i,i) = 1.0;
        A.insert(j,j) = 1.0/mass[i];
        
        B.insert(i,j) = 1.0;

        B.insert(j,i) = -ks[i] - ks[i+1];
        B.insert(j,j) = -cs[i] - cs[i+1];

        if (i>0){
            B.insert(j,i-1) = ks[i];
            B.insert(j,j-1) = cs[i];
        }
        if (i<n-1){
            B.insert(j,i+1) = ks [i+1];
            B.insert(j,j+1) = cs [i+1];
        }
    }

// INTEGRATION
    std::vector<std::vector<double> > coeffs_all_sparse;

    // double t=0;

    //Compute constant zero force, put inside loop if non-constant

    // cout << xv << endl; //DEBUG

    for(int i=0; i<tend/step; i++){
        // std::cout<<"iteration "<<i<<std::endl;
        xv = propagate_euler (xv,step,A,B);
        // xv = xv + step*(A*(B*xv));
        
        // t+=step;

        if((i+1)%freq == 0){
            for (int v=0; v<2*n; v++){// if only wanna track some masses, change here with repr
                Eigen::VectorXd y(npoints);
                for (int j=0;j<npoints;j++){
                    y(j)=xv(v,j);
                }

                // Eigen::VectorXd coe = base_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
                Eigen::VectorXd coe = base_inv*y;
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
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    cout << n<<"-dof spring-mass, dim_unc="<<dim_unc<<", non-intr. sparse, time elapsed : " << time_elapsed << endl << endl;

    // write to file
    std::vector<std::vector<double> > coeffs_all;
    coeffs_all=sparse_to_full(coeffs_all_sparse,idx,dim_unc,degree);

    std::ofstream file;
    file.open ("sm_ni_sparse.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}