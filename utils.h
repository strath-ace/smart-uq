#ifndef UTILS_H
#define UTILS_H

#include <vector>

inline int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
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





#endif // UTILS_H
