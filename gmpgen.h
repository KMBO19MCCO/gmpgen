//
// Created by eszdman on 14.04.23.
//

#ifndef GMPGEN_GMPGEN_H
#define GMPGEN_GMPGEN_H

#include <iostream>
#include <gmpxx.h>
#include <random>
#include <omp.h>
#include <cassert>

using namespace std;

///Computes (a*b - c*d) with precision not worse than 1.5*(unit of the least precision) suggested in Claude-Pierre Jeannerod,
//Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
//Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264
/// \tparam fp_t The data type used for calculations
/// \param a
/// \param b
/// \param c
/// \param d
/// \return calculated result
template<typename fp_t>
fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d);


template<typename fp_t>
class Framework {
public:
    Framework() = default;

    vector<fp_t> roots;
    vector<fp_t> coefficients;

    /// Framework constructor
    /// \param rootsCount number of roots
    explicit Framework(int rootsCount) {
        roots = vector<fp_t>(rootsCount);
        coefficients = vector<fp_t>(rootsCount + 1);
    }


    friend ostream &operator<<(ostream &os, const Framework &framework) {
        os << "size: " << framework.roots.size() << endl;
        os << "seed: " << framework.internalSeed << endl;
        os << "roots: ";
        for (const auto &root: framework.roots) {
            os << root << ' ';
        }
        os << endl;
        os << "coefficients: ";
        for (const auto &coefficient: framework.coefficients) {
            os << coefficient << ' ';
        }
        os << endl;
        return os;
    }
    /// Generating and testing a set of roots according to the specified parameters
    /// \param count number of samples generated
    /// \param rootsCount the number of roots in the polynomial, the number of roots should not exceed 4
    /// \param low minimum of the range of accepted root values
    /// \param high maximum of the range of accepted root values
    /// \param maxDistance maximum distance between roots, the maximum distance should be less than the segment of root generation
    /// \param multipleRoots number of multiple roots
    /// \param testing a function that implements the calculation of roots
    /// \param seed random seed override
    /// \param slow bignum backend
    static void generateBatch(int count, int rootsCount, fp_t low, fp_t high, fp_t maxDistance, int multipleRoots,
                              vector<fp_t> (*testing)(vector<fp_t>), unsigned long long seed = 0, bool slow = false);

private:


    /// Generating roots according to the specified parameters
    /// \param low minimum of the range of accepted root values
    /// \param high maximum of the range of accepted root values
    /// \param maxDistance maximum distance between roots
    /// \param multipleRoots number of multiple roots
    void generate(fp_t low, fp_t high, fp_t maxDistance, int multipleRoots, mt19937_64 &generator);

    /// Calculating the deviation between the original roots and the calculated ones
    /// \param rootsInput the resulting roots
    /// \return a pair of absolute and relative errors
    pair<fp_t, fp_t> deviation(vector<fp_t> rootsInput);

    void generateSlow(fp_t low, fp_t high, fp_t maxDistance, int multipleRoots, mt19937_64 &generator);
};

#endif //GMPGEN_GMPGEN_H
