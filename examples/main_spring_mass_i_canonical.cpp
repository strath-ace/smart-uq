#include "../include/smartuq.h"

// std::vector<Canonical_Polynomial <double> > derivxv (std::vector<Canonical_Polynomial <double> > xv, std::vector <double> ks, std::vector <double> cs, std::vector<double> mass){
//     Eigen::MatrixXd dxv = A*(B*xv);
//     return dxv;
// }std::vector<Canonical_Polynomial <double> >

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

template < class T >
std::vector<Canonical_Polynomial <double> > propagate_euler (std::vector<Canonical_Polynomial <double> > xv, double step, std::vector <double> ks, std::vector <double> cs, T mass){
    
    int n = xv.size()/2;

    std::vector<Canonical_Polynomial <double> > res;
    Canonical_Polynomial<double> term(xv[0].get_nvar(),xv[0].get_degree());

    // //DEBUG
    // clock_t begin = clock();
    
    // x
    for (int i=0; i<n; i++){
        term = xv[i] + step*xv[i+n];
        res.push_back(term);
    }

    // //DEBUG
    // std::cout << "t1=" << ((double)(clock()-begin))/CLOCKS_PER_SEC<<endl;
    // begin=clock();

    // v
    for (int i=0; i<n; i++){
        term = -(ks[i]+ks[i+1])*xv[i] - (cs[i]+cs[i+1])*xv[i+n];
        if (i>0) term += ks[i]*xv [i-1] + cs[i]*xv[i+n-1];
        if (i<n-1) term += ks[i+1]*xv[i+1] + cs[i+1]*xv[i+n+1];
        term *= step;
        term /= mass[i];
        term += xv[i+n];
        res.push_back(term);
    }

    // //DEBUG
    // std::cout << "t2=" << ((double)(clock()-begin))/CLOCKS_PER_SEC<<endl<<endl;

    return res;
}

int main()
{
for (int n=1; n<=12; n++){
// Coded only for uncertainty in states (x and v). It can be in all or just in some of them.

// INPUT

// Algebra
    int degree = 4;

// System
    // n_DOF of the system
    // int n = 20;
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
    // std::vector<double> unc_fs(n,0.0); // NOT IMPLEMENTED    

// Simulation
    double step = 0.005;
    double tend = 50.0;
    int freq = 500; //every how many iterations we save the results

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
    //Creation of vector of ranges
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

    // std::cout<<combination(dim_unc,degree)<<endl;
    //timer
    clock_t begin,end;
    begin=clock();

    // to accelerate multiplication
    Canonical_Polynomial<double>::initialize_M(dim_unc,degree);


    //assign initial status, there will be 2n polynomial representations of dimension dim_unc
    std::vector<Canonical_Polynomial <double> > xv;
    
    // initialize x polynomials
    int count=0; 
    for (int i=0; i<n;i++){
        Canonical_Polynomial<double> poly(dim_unc, degree);
        if (!idx[0].empty() && idx[0][count]==i){
            poly.set_coeffs(count+1,1);
            poly *= (ranges[count][1]-ranges[count][0])/2.0;
            poly += (ranges[count][1]+ranges[count][0])/2.0;
            count++;
        }
        else{
            poly.set_coeffs(0,x0[i]);
        }        
        xv.push_back(poly);
    }
    // initialize v polynomials
    count=0; 
    for (int i=0; i<n;i++){
        Canonical_Polynomial<double> poly(dim_unc, degree);
        if (!idx[1].empty() && idx[1][count]==i){
            poly.set_coeffs(dim[0]+count+1,1);
            poly *= (ranges[dim[0]+count][1]-ranges[dim[0]+count][0])/2.0;
            poly += (ranges[dim[0]+count][1]+ranges[dim[0]+count][0])/2.0;
            count++;
        }
        else{
            poly.set_coeffs(0,v0[i]);
        }        
        xv.push_back(poly);
    }
    // initialize m polynomials
    std::vector<Canonical_Polynomial <double> > mass_polys;
    bool massisunc = (dim[2]!=0);
    if (massisunc){
        count=0; 
        for (int i=0; i<n;i++){
            Canonical_Polynomial<double> poly(dim_unc, degree);
            if (!idx[2].empty() && idx[2][count]==i){
                poly.set_coeffs(dim[0]+dim[1]+count+1,1);
                poly *= (ranges[dim[0]+dim[1]+count][1]-ranges[dim[0]+dim[1]+count][0])/2.0;
                poly += (ranges[dim[0]+dim[1]+count][1]+ranges[dim[0]+dim[1]+count][0])/2.0;
                count++;
            }
            else{
                poly.set_coeffs(0,mass[i]);
            }        
            mass_polys.push_back(poly);
        }
    }

    // HERE WE USUALLY  BUILD SYSTEM MATRICES BUT IT SEEMS CUMBERSOME TO MAKE EIGEN ACCEPT POLYNOMIAL CLASS.
    // INSTEAD WE MODIFY THE DERIVXV AND PROPAGATE FUNCTIONS TO AVOID MATRICES AND WORK WITH K, C, M VECTORS

// INTEGRATION
    std::vector<std::vector<double> > coeffs_all;
    // double t=0;
    //Compute constant zero force, put inside loop if non-constant NOT IMPLEMENTED

    for(int i=0; i<tend/step; i++){
        // if (i%100==0) std::cout<<"iteration "<<i<<std::endl;

        if (massisunc) xv = propagate_euler (xv,step,ks,cs,mass_polys);
        else xv = propagate_euler (xv,step,ks,cs,mass);

        // t+=step;

        if((i+1)%freq == 0){
            for (int v=0; v<2*n; v++){// if only wanna track some masses, change here with repr
                std::vector<double> coeffs = xv[v].get_coeffs();
                coeffs_all.push_back(coeffs);
            }
        }
    }
    //timer

    //deallocate M for next thing
    Canonical_Polynomial<double>::delete_M();

    end=clock();
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    cout << n<<"-dof spring-mass, dim_unc="<<dim_unc<<", intr.canonical, time elapsed : " << time_elapsed << endl << endl;
    
    // write to file
    std::ofstream file;
    file.open (("sm_i_canonical_n"+patch::to_string( n)+".out").data());
    // file.open ("sm_i_canonical.out);
    for(int k=0; k<coeffs_all.size(); k++){
        for(int kk=0; kk<coeffs_all[k].size(); kk++){
            file  << setprecision(16) << coeffs_all[k][kk] << " ";
        }
        file << "\n";
    }
    file.close();
}
}

template class std::vector< Canonical_Polynomial <double> >;
template class std::vector< double >;

