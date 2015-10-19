#include "main_list.h"

void main_tests(){

//    //TEST FUNCTIONS UNIVARIATE
//    int nvar = 1;
//    for(int degree = 2; degree <11 ; degree++){
//    std::vector<Chebyshev_Polynomial<double> > x,f;
//    for(int i=0; i<8; i++){
//        x.push_back(Chebyshev_Polynomial<double>(nvar,degree));
//        x[i].set_coeffs(1,1);
//    }

//    std::ofstream file;
//    std::string filename = "approx_"+patch::to_string(degree)+".out";
//    file.open (filename.data());

//    x[0] = (4.0-3.0)/2.0*x[0] + (4.0+3.0)/2.0;
//    Chebyshev_Polynomial<double> f1 = sin(x[0]);
//    f.push_back(f1);

//    x[1] = (1.0)/2.0*x[1] + (1.0)/2.0;
//    Chebyshev_Polynomial<double> f2 = exp(1.0/cos(x[1]));
//    f.push_back(f2);

//    x[2] = (1.0)/2.0*x[2] + (1.0)/2.0;
//    Chebyshev_Polynomial<double> f3 = exp(x[2])/(log(2.0+x[2])*cos(x[2]));
//    f.push_back(f3);

//    Chebyshev_Polynomial<double> f4 = sin(exp(x[3]));
//    f.push_back(f4);

//    x[4] = (1.0)/2.0*x[4] + (-1.0)/2.0;
//    Chebyshev_Polynomial<double> f5 = sqrt(x[4]+1.0001);
//    f.push_back(f5);

//    x[5] = (1.0)/2.0*x[5] + (-1.0)/2.0;
//    Chebyshev_Polynomial<double> f6 = sqrt(x[5]+1.0001)*sin(x[5]);
//    f.push_back(f6);

//    Chebyshev_Polynomial<double> f7 = 1.0/(1.0+4.0*x[6]*x[6]);
//    f.push_back(f7);

//    Chebyshev_Polynomial<double> f8 = pow(sin(x[7]),2) + pow(cos(x[7]),2);
//    f.push_back(f8);

//    for(int i=0; i<8;i++){
//        std::vector<double> coeffs = f[i].get_coeffs();
//        for(int k=0; k<coeffs.size(); k++){
//            file << setprecision(16)<<coeffs[k] << " ";
//        }
//        file << "\n";
//    }
//    }

        //TEST FUNCTIONS MULTIVARIATE
        int nvar = 2;
        for(int degree = 5; degree <=5 ; degree++){
        std::vector<Chebyshev_Polynomial<double> > x, f;
        for(int i=0; i<nvar; i++){
            x.push_back(Chebyshev_Polynomial<double>(nvar,degree));
            x[i].set_coeffs(1+i,1);
        }

        // std::ofstream file;
        // std::string filename = "approx_mv_"+patch::to_string(degree)+".out";
        // file.open (filename.data());

        x[0] = (4.0-3.0)/2.0*x[0] + (4.0+3.0)/2.0;    // x in [3.0,4.0]
        //x[1] = (1.0-0.0)/2.0*x[1] + (1.0+0.0)/2.0;    // y in [0,1]
        //x[2] = (2.0-1.0)/2.0*x[2] + (2.0+1.0)/2.0;    // z in [1,2]
        Chebyshev_Polynomial<double> f1 = 2.0*x[0]*x[0]*x[0];//x[0]*x[0]*x[0]+5.0*x[0]*x[1]+x[2]*x[2]*x[1]+1.0;
        //Chebyshev_Polynomial<double> f2 = 5.0*x[0]*x[1];//pow(x[2],3)+2.0*x[0]*pow(x[2],2)+pow(x[1],2);
        cout << 1/f1 << endl;
        //f.push_back(f1);
        //f.push_back(f2);

        Chebyshev_Polynomial<double> ff(nvar,degree);
//        ff = f1/f2;
//        f.push_back(ff);

        ff = 1.0;//sin(f1)*cos(f2);
        f.push_back(ff);

//        ff = f1*f2;
//        f.push_back(ff);

//        ff = sqrt(f1+f2);
//        f.push_back(ff);

        ff = sin(f1);
        f.push_back(ff);

        ff = cos(f1);
        f.push_back(ff);

        // for(int i=0; i<f.size();i++){
        //     std::vector<double> coeffs = f[i].get_coeffs();
        //     for(int k=0; k<coeffs.size(); k++){
        //         file << setprecision(16)<<coeffs[k] << " ";
        //     }
        //     file << "\n \n";
        // }

        }


//    //TEST INTEGRATION

//    int nvar = 1;
//    std::ofstream file;
//    int degree = 10;
//    file.open ("approx_ode_10.out");

//    for(int tend = 1; tend<9; tend++){

//        std::vector<Chebyshev_Polynomial<double> >  x;
//        x.push_back(Chebyshev_Polynomial<double>(nvar,degree));
//        x[0].set_coeffs(1,1);

//        double step = 0.01;

//        std::vector<double> ranges(2);
//        ranges[0] = 1.0;
//        ranges[1] = 2.0;
//        //translation  [-1,1] ----> [a,b]
//        x[0] = (ranges[1]-ranges[0])/2.0*x[0] + (ranges[1]+ranges[0])/2.0;

//        //assign initial status

//        //perform integration
//        for(int i=0; i<tend/step; i++){
//            std::cout<<"iteration "<<i<<std::endl;

//            x = rk4<double>(f,x,step);

//        }

//        std::vector<double> coeffs = x[0].get_coeffs();
//        for(int k=0; k<coeffs.size(); k++){
//            file << coeffs[k] << " ";
//        }
//        file << "\n";
//    }

//    file.close();
}
