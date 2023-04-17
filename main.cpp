//
// Created by eszdman on 14.04.23.
//
#include "gmpgen.h"

template<typename fp_t>
vector<fp_t> solve(std::vector<fp_t> coefficients) {
    vector<fp_t> output(coefficients.size() - 1);
    fp_t A, B, C;
    //Coefficients should be in ABC order
    A = coefficients[0];
    B = coefficients[1];
    C = coefficients[2];
    fp_t d = pr_product_difference(B,B, 4*A, C);
    d = std::max(d, static_cast<fp_t>(numeric_limits<fp_t>::epsilon()));
    if (d > numeric_limits<fp_t>::epsilon()) {
        output[0] = (-B + sqrt(d)) / (2 * A);
        output[1] = (-B - sqrt(d)) / (2 * A);
    } else if (d == numeric_limits<fp_t>::epsilon()) {
        output[0] = (-B / 2 * A);
        output[1] = (-B / 2 * A);
    }
    return output;
}

int main() {
    Framework<float>::generateBatch(1, 2, -1, 1, 1e-4, 2, solve, 0, false);
    return 0;
}
