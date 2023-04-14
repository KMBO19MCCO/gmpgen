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



class Framework {
public:
    vector<mpf_class> roots;
    vector<float> rootsReal;
    vector<mpf_class> coefficients;
    vector<float> coefficientsReal;
    unsigned long long internalSeed;

    explicit Framework(int rootsCount, unsigned long long seed = 0) {
        roots = vector<mpf_class>(rootsCount);
        rootsReal = vector<float>(rootsCount);
        coefficients = vector<mpf_class>(rootsCount + 1);
        coefficientsReal = vector<float>(rootsCount + 1);

        internalSeed = randomDevice();
        if (seed) internalSeed = seed;
        generator = mt19937_64(internalSeed);
        cout.precision(50);
    }

    void generate(float low, float high, float dist);

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

    static void generateBatch(int count, int rootsCount, float low, float high, float maxDistance, vector<float> (*testing)(vector<float>));

private:
    random_device randomDevice;
    mt19937_64 generator;

    pair<float,float> deviation(vector<float> roots);
};

#endif //GMPGEN_GMPGEN_H
