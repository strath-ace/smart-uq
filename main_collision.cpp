#include "main_list.h"

void main_collision_avoidance(){
    std::ofstream file;
    file.open ("results2.out");

    double pi = 3.141592653589793;

    //algebra params
    int degree = 4;
    int nvar = 4;
    //integration params
    double step = 0.01;
    double tend = 1.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    std::vector<double> x(4), unc(4);
        x[0] = 2.362068277608905;
        x[1] = -2.662198239173287;
        x[2] = -0.169213456203144;
        x[3] = 0.215089310029571;

//    x[0] = 1.0;
//    x[1] = -5.0;
//    x[2] = 0.0;
//    x[3] = 1.0;

    unc[0] = 0.01;
    unc[1] = 0.01;
    unc[2] = 0.01;
    unc[3] = 0.01;

    //2BP pictures
    ranges[0][0] = x[0]-unc[0];
    ranges[0][1] = x[0]+unc[0];

    ranges[1][0] = x[1]-unc[1];
    ranges[1][1] = x[1]+unc[1];

    ranges[2][0] = x[2]-unc[2];
    ranges[2][1] = x[2]+unc[2];

    ranges[3][0] = x[3]-unc[3];
    ranges[3][1] = x[3]+unc[3];

    //

    std::vector<Chebyshev_Polynomial<double> > x0,param0;
    for(int i=0; i<nvar; i++){
        x0.push_back(Chebyshev_Polynomial<double>(nvar,degree));
        x0[i].set_coeffs(i+1,1);
    }

    std::vector<Chebyshev_Polynomial<double> > res;
    for(int i=0; i<nvar; i++){
        res.push_back(Chebyshev_Polynomial<double>(nvar,degree));
    }

    //translation  [-1,1] ----> [a,b]
    for(int i=0; i<nvar; i++){
        x0[i] = (ranges[i][1]-ranges[i][0])/2.0*x0[i] + (ranges[i][1]+ranges[i][0])/2.0;
    }

    //assign initial status
    for(int i=0; i<nvar; i++){
        res[i] = x0[i];
    }

    //perform integration
    for(int i=0; i<tend/step; i++){
        std::cout<<"iteration "<<i<<std::endl;

        res = rk4<double>(f,res,param0,step);

        for(int j=0; j<nvar; j++){
            std::vector<double> coeffs = res[j].get_coeffs();
            for(int k=0; k<res[j].get_coeffs().size(); k++){
                file << coeffs[k] << " ";
            }
            file << "\n";
        }
    }
    file.close();
}

