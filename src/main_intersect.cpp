#include "main_list.h"
#include <sstream>

void main_multiphase(){
    
    //algebra params
    int degree = 4;
    int nvar = 4;
    int nphases = 3; // number of intersection calculations
    int nprop = 3; // number of parallel propagations
    	
    //integration params
    double step = 0.01;
    double sma = 2;
    double tend = 2.0*M_PI/pow(sma,-3.0/2.0);
    double e = 0.5;

    std::vector<double> x(nvar);
    std::vector<double> unc_x(nvar);

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = sqrt(1+e);

    unc_x[0] = 0.0;
    unc_x[1] = 0.0;
    unc_x[2] = 0.005;
    unc_x[3] = 0.005;

    // initialize randomly nprop values for the velocity
    std::vector<std::vector<double> > states;
    std::vector<std::vector<double> > unc_states;
    for(int i=0; i<nprop; i++){
        states.push_back(std::vector<double>(nvar));
        states[i][0] = x[0];
        states[i][1] = x[1];
    	for(int j=2; j<4; j++){
		double vel = (double)rand() / RAND_MAX;
        	states[i][j] = x[j]-unc_x[j] + vel * (2*unc_x[j]);
    	}

	unc_states.push_back(std::vector<double>(nvar));
        unc_states[i][0] = unc_x[0];
        unc_states[i][1] = unc_x[1];
	unc_states[i][2] = unc_x[2]/nprop;
	unc_states[i][2] = unc_x[3]/nprop;
    }

    
    for(int n=0; n<nprop; n++){
	//LOOP ON PROPAGATIONS
	std::ofstream file;
	std::ostringstream stm ;
        stm << n ;
	std::string name_file = "intersection_" + stm.str() + ".out";
	file.open (name_file.data());
	
	//setting ranges expansion
    	std::vector<std::vector<double> > ranges_states;
    	for(int i=0; i<nvar; i++){
        	ranges_states.push_back(std::vector<double>(2));
    	}
    	for(int i=0; i<nvar; i++){
        	ranges_states[i][0] = states[n][i]-unc_states[n][i];
        	ranges_states[i][1] = states[n][i]+unc_states[n][i];
    	}
	//initialize initial state in the tchnebycheff base
	std::vector<Chebyshev_Polynomial<double> > x0;
        for(int i=0; i<nvar; i++){
        	x0.push_back(Chebyshev_Polynomial<double>(nvar,degree));
        	x0[i].set_coeffs(i+1,1);
    	}
        //allocate memory for final state vector of polynomials
	std::vector<Chebyshev_Polynomial<double> > res;	
        for(int i=0; i<nvar; i++){
        	res.push_back(Chebyshev_Polynomial<double>(nvar,degree));
    	}
        //translation  [-1,1] ----> [a,b]
        for(int i=0; i<nvar; i++){
        	x0[i] = (ranges_states[i][1]-ranges_states[i][0])/2.0*x0[i] + (ranges_states[i][1]+ranges_states[i][0])/2.0;
        }
    	//assign initial status
    	for(int i=0; i<nvar; i++){
        	res[i] = x0[i];
    	}

	std::vector<std::vector<double> > coeffs_all;
    	try{
            //perform integration
            for(int i=0; i<tend/step; i++){
            	std::cout<<"iteration "<<i<<std::endl;
		std::vector<Chebyshev_Polynomial<double> > param;
		param.push_back(Chebyshev_Polynomial<double>(1,1));

            	res = rk4<double>(f,res,param,step);

            	//save values to be printed
            	if((i+1)%100 == 0){
                	for(int j=0; j<nvar; j++){
                    		std::vector<double> coeffs = res[j].get_coeffs();
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
        	}
        	file.close();
        	exit(EXIT_FAILURE);
    	}


    	for(int k=0; k<coeffs_all.size(); k++){
        	for(int kk=0; kk<coeffs_all[k].size(); kk++)
            		file << setprecision(16) << coeffs_all[k][kk] << " ";
        	file << "\n";
    	}
    	file.close();


    }//close loop on multiple propagations


}


