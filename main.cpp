//
// Created by eszdman on 14.04.23.
//
#include "gmpgen.h"
template<typename fp_t>
fp_t dscrmt(fp_t A,fp_t B,fp_t C){
    return B*B - 4*A*C;
}

template<typename fp_t>
vector<float> solve(std::vector<fp_t> coefficients){
    vector<fp_t> output(coefficients.size()-1);
    fp_t A,B,C;
    //Coefficients should be in ABC order
    A = coefficients[0];
    B = coefficients[1];
    C = coefficients[2];
    fp_t d = dscrmt(A,B,C);
    d = std::max(d,static_cast<fp_t>(std::numeric_limits<fp_t>::epsilon()));
    output[0] = (-B - sqrt(d))/(2*A);
    output[1] = (-B + sqrt(d))/(2*A);
    return output;
}



int main() {
    Framework::generateBatch(1000'0000, 2, -1., 1, 1e-1, solve<float>);
    return 0;
}
