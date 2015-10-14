#include <iostream>     // std::cout
#include <algorithm>    // std::next_permutation, std::sort
#include <vector>       // std::vector
#include <fstream>      // I/O on file
#include <sstream>
#include <string.h>
#include <iterator>
#include <math.h>
#include <stdio.h>

using namespace std;

vector<int> sumVector( vector<int> a, vector<int> b);
vector<int> array2vector(int arr[], int n);
void printArray(int arr[], int size);
void printVector(vector<int>& v, int level, ofstream& myfile);
void print_rows(vector< vector<int> > vv, string filnam);
vector< vector<int> > part_recursive(int n, vector<int>& v, int level, vector< vector<int> > vv, int dim);
vector< vector<int> > partitions(int N, int dim);
vector< vector<int> > permutations( vector<int> v, int dim, vector< vector<int> > vv );
int emme(int i);
vector< vector<int> > smolyak_recursive(int n, int dim, int lower[], int upper[],vector<int> v, vector< vector<int> > vv);
vector< vector<int> > smolyak_indices(vector<int> v, int dim, vector< vector<int> > vv);
vector< vector<int> > sparse_grid_index( int dim, int mu);
vector< vector<int> > initialize_J( int dim, int deg);
vector< vector<int> > initialize_N( int dim, int deg, vector< vector<int> > Jtab);
int vector_to_idx(vector<int> k, int dim, vector< vector<int> > Ntab);
vector<int> idx_to_vector( int idx, int deg, int dim, vector< vector<int> > Ntab);
