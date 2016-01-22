#include "../include/smartuq.h"

Eigen::MatrixXd derivxv(Eigen::MatrixXd xv, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
    Eigen::MatrixXd dxv = A*(B*xv);
    return dxv;
}

Eigen::MatrixXd propagate_euler (Eigen::MatrixXd xv, double step, Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B, std::vector < std::vector <int> > idx, std::vector <int> dim, std::vector<std::vector<double> > par){
    Eigen::MatrixXd xv_next(xv.rows(),xv.cols());
    if (dim[2]==0){
        xv_next = xv + step*derivxv(xv,A,B);
    }
    else{
        for (int i=0; i<xv.cols(); i++){
            // assign nominal A
            Eigen::SparseMatrix<double> Ai = A;
            // modify A
            for (int j=0; j<dim[2]; j++){
                int jj=idx[2][j];
                Ai.coeffRef(jj,jj)=1.0/par[i][dim[0]+dim[1]+j];
            } 
            // propagate only this point   
            xv_next.col(i) = xv.col(i) + step*derivxv(xv.col(i),Ai,B);
        }
    }
    return xv_next;
}

int main()
{

// Coded only for uncertainty in states (x and v). It can be in all or just in some of them.

// INPUT

// Algebra
    int degree = 4;

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
    // std::vector<double> fs(n,0.0); // amplitudes of the exciting forces NOT IMPLEMENTED

    // // Uncertainty in parameters 
    std::vector<double> unc_mass(n,0.05);
    // std::vector<double> unc_ks(n+1,0.0); // NOT IMPLEMENTED
    // std::vector<double> unc_cs(n+1,0.0); // NOT IMPLEMENTED
    // std::vector<double> unc_fs(n,0.0);   // NOT IMPLEMENTED   

// Simulation
    double step = 0.005;
    double tend = 50.0;
    int freq = 2500; //every how many iterations we save the results

// INITIALISATIONS

    std::vector < std::vector <int> > idx(3);
    std::vector <int> dim(3,0);
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

        if (fabs(unc_mass[i])>ZERO) {
            idx[2].push_back(i);
            dim[2]++;
            dim_unc++;
        }

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

// INITIALISATIONS
    // Number of points and coefficients
    int npoints=combination(dim_unc,degree);
    int ncoeffs=combination(dim_unc,degree);
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

    for (int i=0; i<dim[2];i++){
        std::vector<double> r(2);
        r[0] = mass[idx[2][i]]-unc_mass[idx[2][i]];
        r[1] = mass[idx[2][i]]+unc_mass[idx[2][i]];
        ranges.push_back(r);
    }


    //timer
    clock_t begin,end;
    begin=clock();

    //sampling and evaluating basis, there will be 2n polynomial representations of dimension dim_unc
    sampling::lhs<double> lhs_gen(dim_unc,npoints);
    Eigen::MatrixXd base_matrix (npoints,ncoeffs);

    std::vector<std::vector<double> > par;

    for(int i=0; i<npoints; i++){
        std::vector<double> par_aux=lhs_gen();
        std::vector<double> par_next;
        // translation to real range for x0 and to -1,1 for evaluation of basis
        for (int j=0; j<dim_unc; j++){
            par_next.push_back(ranges[j][0]+par_aux[j]*(ranges[j][1]-ranges[j][0]));
            par_aux[j]*=2;
            par_aux[j]-=1;
        }

        par.push_back(par_next);

        //evaluation of the base
        Chebyshev_Polynomial<double> base_poly(dim_unc,degree,1.0);
        base_matrix(i,0)=1.0;
        for (int j=1;j<ncoeffs;j++){
            base_poly.set_coeffs(j-1,0.0);
            base_poly.set_coeffs(j,1.0);
            base_matrix(i,j)=base_poly.evaluate(par_aux);
        }
    }

    // //DEBUG
    // for (int i=0;i<par.size();i++){
    //     for (int j=0;j<par[i].size();j++){
    //         std::cout<< par[i][j]<<"    ";
    //     }
    //     std::cout<<endl<<endl;
    // }

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
    std::vector<std::vector<double> > coeffs_all;

    // double t=0;

    //Compute constant zero force, put inside loop if non-constant NOT IMPLEMENTED

    // cout << xv << endl; //DEBUG

    for(int i=0; i<tend/step; i++){
        // if(i%100==0) std::cout<<"iteration "<<i<<std::endl;

        xv = propagate_euler (xv,step,A,B,idx,dim,par);
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
                coeffs_all.push_back(coeffs);
            }
        }
    }
    //timer
    end=clock();
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    cout << n<<"-dof spring-mass, dim_unc="<<dim_unc<<", non-intr. full, time elapsed : " << time_elapsed << endl << endl;
    
    // write to file
    std::ofstream file;
    file.open ("sm_ni.out");
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}




