//
// Created by eszdman on 14.04.23.
//

#ifndef GMPGEN_GMPGEN_H
#define GMPGEN_GMPGEN_H

#include <iostream>
#include <gmpxx.h>
#include <random>

using namespace std;

/* computes (a*b - c*d) with precision not worse than 1.5*(unit of the least precision) suggested in Claude-Pierre Jeannerod,
Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264 */
template<typename fp_t>
fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d);


template<typename fp_t>
class Framework {
public:
    vector<mpf_class> roots;
    vector<fp_t> rootsReal;
    vector<mpf_class> coefficients;
    vector<fp_t> coefficientsReal;
    unsigned long long internalSeed;

    explicit Framework(int rootsCount, unsigned long long seed = 0) {
        roots = vector<mpf_class>(rootsCount);
        rootsReal = vector<fp_t>(rootsCount);
        coefficients = vector<mpf_class>(rootsCount + 1);
        coefficientsReal = vector<fp_t>(rootsCount + 1);

        internalSeed = randomDevice();
        if (seed) internalSeed = seed;
        generator = mt19937_64(internalSeed);
        cout.precision(50);
    }

    void generate(fp_t low, fp_t high, fp_t dist);

    friend ostream &operator<<(ostream &os, const Framework &framework) {
        //framework.generator.seed();
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

    static void generateBatch(int count, int rootsCount, fp_t low, fp_t high, fp_t maxDistance,
                              vector<fp_t> (*testing)(vector<fp_t>));

private:
    random_device randomDevice;
    mt19937_64 generator;

    pair<fp_t, fp_t> deviation(vector<fp_t> rootsInput);
};

#endif //GMPGEN_GMPGEN_H
