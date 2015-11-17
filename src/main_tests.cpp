#include "main_list.h"
#include "utils.h"

void main_tests(){

    // int nvar, degree;

    // cout << "degree max =  ";
    // cin >> degree;
    // cout << endl;

    // cout << "number of variables =  ";
    // cin >> nvar;
    // cout << endl;

    // std::ofstream file;
    // std::string filename = "canonical_poly_"+patch::to_string(degree)+"_"+patch::to_string(nvar)+".m";
    // file.open (filename.data());
    // file << "function y = canonical_poly(x)" << "\n" << "% Canonical Polynomial basis of up to order "<<degree<<" in "<<nvar<<" variables" << "\n" << "\n" << "y = [ 1,..." << "\n";

    // Canonical_Polynomial<double> p(nvar,degree);
    // for (int deg=1;deg<=degree;deg++){
    //     for(int i=0; i<p.get_J()[nvar][deg]; i++){
    //         std::vector<int> row = p.get_row(i,deg);
    //         for (int j=0;j<nvar;j++){
    //             if (row[j]==0) file << 1;
    //             else if (row[j]==1) file << "x(" << j+1 << ")";
    //             else file << "x(" << j+1 << ").^"<< row[j];
    //             if (j!= nvar-1) file << " .* ";
    //         }
    //         if (deg==degree && i==p.get_J()[nvar][deg]-1) file << "  ];"<< "\n" << "end";
    //         else file << ",..." << "\n";
    //     }
    // }
    // file.close();
} 