#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>


template <class T>
T inverse(T x){
    return 1.0/x;
}


//MATH STUFFS
inline int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline int combination(int n, int k)
{
    int max = std::max(n,k);
    int min = std::min(n,k);
    int res = 1;
    int j = 1;

    while(j<=min){
        res *= max+j;
        j++;
    }

    return res/factorial(min);
}


inline void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, int count){
    if (count < item.size()){
        for (int i = 0; i < values.size(); i++) {
            item[count] = values[i];
            int tmp_count = count + 1;
            rep(res, values, item, tmp_count);
        }
    }else{
        res.push_back(item);
    }
}


inline void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res){
    res.clear();

    std::vector<int> item(k);
    rep(res, values, item, 0);
}



//LAPACK METHOD
//    std::vector<T> coeffs = other.get_coeffs();
//    int n = coeffs.size();
//    int nrhs = 1;
//    double B[n][n];
//    double b[1][n];
//    int lda = n;
//    int ldb = n;
//    int ipiv[n];
//    int info;

//    for(int i=0; i<n; i++){
//        b[0][i] = 0.0;
//        for(int j=0; j<=i; j++){
//            if(j==0){//first row and first diagonal
//                B[i][0] = coeffs[i];
//                B[j][i] =  B[i][j];
//            }
//            else if (i==j){//diagonal
//                if(2*i < n)
//                    B[i][i] = 2.0*coeffs[0]+coeffs[2*i];
//                else
//                    B[i][i] = 2.0*coeffs[0];
//            }
//            else{//lower diagonal
//                if((i+j)<n)
//                    B[i][j] = coeffs[fabs(i-j)]+coeffs[i+j];
//                else
//                    B[i][j] = coeffs[fabs(i-j)];
//                B[j][i] =  B[i][j];
//            }
//        }
//    }

//    b[0][0] = 2.0;

//    dgesv_(&n, &nrhs, &B[0][0], &lda, ipiv, &b[0][0], &ldb, &info);

//    // Check for success
//    if(info == 0)
//    {
//       std::vector<T> res_coeffs(n);
//       for(int i=0; i<n; i++)
//           res_coeffs[i] = b[0][i];
//       res_coeffs[0] /= 2.0;
//       res.set_coeffs(res_coeffs);
//    }
//    else
//    {
//       // Write an error message
//       std::cout << "LAPACK dgesv returned error " << info << "\n";
//       exit(EXIT_FAILURE);
//    }


#endif // UTILS_H
