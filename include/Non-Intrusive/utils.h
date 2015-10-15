# include <string>
# include <ctime>
# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
#include <vector>       // std::vector
#include <iterator>
#include <string.h>

using namespace std;

int i4_min ( int i1, int i2 );
int i4_max ( int i1, int i2 );
string i4_to_string ( int i4, string format );
double r8_epsilon ( );
void timestamp ( );
void print_double_vector(vector<double> v, string filnam, string header);
void print_double_matrix(vector< vector<double> > vv, string filnam, string header);
