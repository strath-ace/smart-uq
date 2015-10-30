#include "main_list.h"
#include "utils.h"

void main_tests(){

    srand(time(NULL));

    int nvar=1;
    int deg=4;
    Newton_Polynomial<double> p(nvar,deg);

    cout << p << endl;
    std::vector<double> nodes(5), coeffs(5);
    nodes[0]=-1.5;
    nodes[1]=-.75;
    nodes[2]=0;
    nodes[3]=.75;
    nodes[4]=1.5;
    coeffs[0]=-14.1014;
    coeffs[1]=17.5597;
    coeffs[2]=-10.8784;
    coeffs[3]=4.83484;
    coeffs[4]=0;

    p.set_nodes(nodes);
    p.set_coeffs(coeffs);
    
    for (int i=0; i<nodes.size(); i++) cout << nodes[i] << "    ";
    cout << endl;
    
    cout << p.evaluate(0)<<endl;
    cout << p.evaluate(-.75)<<endl;
    cout << p.evaluate(1.5)<<endl;

}