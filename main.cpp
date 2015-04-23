#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "chebyshev_polynomial.h"
#include "elementary_functions.h"
#include "integrators.h"
#include "f.h"


#include <iterator>
#include <algorithm>


using namespace std;

void main_collision_avoidance(){
    std::ofstream file;
    file.open ("results2.out");

    double pi = 3.141592653589793;

    //algebra params
    int degree = 4;
    int nvar = 4;
    //integration params
    double step = 0.01;
    double tend = 5.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    std::vector<double> x(4), unc(4);
        x[0] =  2.362091431708905;
        x[1] =  -2.662050550712367;
        x[2] =  -0.169170090722185;
        x[3] = 0.215073962962989;
//    x[0] = 1.0;
//    x[1] = -5.0;
//    x[2] = 0.0;
//    x[3] = 1.0;

    unc[0] = 0.01;
    unc[1] = 0.5;
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

    std::vector<Chebyshev_Polynomial<double> > x0;
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

    double lb = ranges[3][0], ub = ranges[3][1];

    std::vector<Chebyshev_Polynomial<double> > f_composed;
    f_composed = f(x0);

    //perform integration
    for(int i=0; i<tend/step; i++){
        std::cout<<"iteration "<<i<<std::endl;
        std::vector<Chebyshev_Polynomial<double> > tmp_res;
        for(int j=0; j<nvar; j++){
            tmp_res.push_back(res[j]);
        }
//        std::vector<Chebyshev_Polynomial<double> > xstart;
//        for(int i=0; i<nvar; i++){
//            xstart.push_back(Chebyshev_Polynomial<double>(nvar,degree));
//            xstart[i].set_coeffs(i+1,1);
//        }
//        //translation  [-1,1] ----> [a,b]
//        for(int i=0; i<nvar; i++){
//            xstart[i] = (ranges[i][1]-ranges[i][0])/2.0*xstart[i] + (ranges[i][1]+ranges[i][0])/2.0;
//        }

        //res = euler<double>(f,res,step);
        //////////////////////////////////////
        for(int j=0;j<nvar;j++)
            res[j] = tmp_res[j] + step*f_composed[j].composition(tmp_res);
        //////////////////////////////////////


//        bool scale_range = false;
//        double delta = 0.0;
//        for(int j=0; j<nvar; j++){
//            double r = tmp_res[j].get_range();
//            if(r > ub && (j==6)){
//                scale_range = true;
//                delta = r - ub;
//                ub = ub + delta;
//                lb = lb - delta;
//                xstart[j] = Chebyshev_Polynomial<double>(nvar,degree);
//                xstart[j].set_coeffs(j+1,1);
//                xstart[j] = (ub-lb)/2.0*xstart[j] + (ub+lb)/2.0;
//            }
//        }
//        if(scale_range){
//            std::vector<Chebyshev_Polynomial<double> > f_scaled;
//            f_scaled = f(xstart);

//            for(int j=0;j<nvar;j++)
//                res[j] = tmp_res[j] + step*f_scaled[j].composition(tmp_res);
//        }

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


void main_2BP(){
    std::ofstream file;
    file.open ("results.out");

    double pi = 3.141592653589793;

    //algebra params
    int degree = 4;
    int nvar = 4;
    //integration params
    double step = 0.01;
    double tend = 5.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    //2BP pictures
    ranges[0][0] = 1.0;
    ranges[0][1] = 1.1;

    ranges[1][0] = -1.0;
    ranges[1][1] = 1.0;

    ranges[2][0] = 0.0;
    ranges[2][1] = 0.1;

    ranges[3][0] = sqrt(1/pow(ranges[0][1],3));
    ranges[3][1] = sqrt(1/pow(ranges[0][0],3));

    std::vector<Chebyshev_Polynomial<double> > x0;
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
        res = euler<double>(f,res,step);

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

void main_vanderpol()
{
    std::ofstream file;
    file.open ("results.out");

    //algebra params
    int degree = 10;
    int nvar = 2;
    //integration params
    double step = 0.01;
    double tend = 5.0;

    std::vector<std::vector<double> > ranges;
    for(int i=0; i<nvar; i++){
        ranges.push_back(std::vector<double>(2));
        ranges[i][0] = -1.0; ranges[i][1] = 1.0;
    }

    ranges[0][0] = -2.0; ranges[0][1] = 2.0;
    ranges[1][0] = -2.0; ranges[1][1] = 2.0;

    std::vector<Chebyshev_Polynomial<double> > x0;
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
        std::vector<Chebyshev_Polynomial<double> > tmp_res;
        for(int j=0; j<nvar; j++){
            tmp_res.push_back(res[j]);
        }

        res = euler<double>(f,res,step);

        bool scale_range = false;
        for(int j=0; j<nvar; j++){
            double r = res[j].get_range();
            if(r > ranges[j][1]){
                scale_range = true;
                x0[j] = Chebyshev_Polynomial<double>(nvar,degree);
                x0[j].set_coeffs(j+1,1);
                x0[j] = r*x0[j];
            }
        }
        if(scale_range){
            std::vector<Chebyshev_Polynomial<double> > f_scaled;
            f_scaled = f(x0);

            for(int j=0;j<nvar;j++)
                res[j] = tmp_res[j] + step*f_scaled[j].composition(tmp_res);
        }


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

int main_orbit_determination(){
    Chebyshev_Polynomial<double> x(1,10);
    x.set_coeffs(1,1);
    Chebyshev_Polynomial<double> y1 = (4.0-x)*(4.0-x)*(5.0+x);
    Chebyshev_Polynomial<double> y2 = x*x-1.0;

    Chebyshev_Polynomial<double> z = y1*y1+y2*y2-1;

    std::cout<<z;

}

int main_test_composition(){
    Chebyshev_Polynomial<double> x(1,10), y(1,10);
    x.set_coeffs(1,1);
    y.set_coeffs(1,1);
    Chebyshev_Polynomial<double> f = x*x + 1.0;
    Chebyshev_Polynomial<double> g = 2.0*y - 3.0;
    std::vector<Chebyshev_Polynomial<double> > gg;
    gg.push_back(g);

    Chebyshev_Polynomial<double> res = 4.0*y*y - 12.0*y + 10.0;
    std::cout << f.composition(gg);
    std::cout << res;
}

int main()
{
    cout << "Welcome to Chebyshev Algebra!" << endl;

    //    test divisione
    //    Chebyshev_Polynomial<double> x(1,50);
    //    x.set_coeffs(1,1);
    //    Chebyshev_Polynomial<double> f = (4.0-x)*(4.0-x)*(5.0+x);
    //    Chebyshev_Polynomial<double> g = 1.0/f;
    //    std::cout<<g<<std::endl;


    //main_2BP();
    //main_collision_avoidance();
    //main_vanderpol();
    //main_orbit_determination();
    main_test_composition();

    return 0;
}
