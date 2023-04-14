#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include <omp.h>
#include "gmpgen.h"

using namespace std;

/* computes (a*b - c*d) with precision not worse than 1.5*(unit of the least precision) suggested in Claude-Pierre Jeannerod,
Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264 */
template<typename fp_t>
fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d) {
    auto tmp = d * c;
    return fma(a, b, -tmp) + fma(-d, c, tmp);
}

pair<float, float> Framework::deviation(vector<float> rootsInput) {
    vector<float> tempRoots = vector<float>(roots.size());
    vector<float> tempRootsInput = vector<float>(roots.size());
    for (auto i = 0; i < roots.size(); ++i) {
        tempRoots[i] = static_cast<float>(roots[i].get_d());
        tempRootsInput[i] = rootsInput[i];
    }
    float maxDev = 0.f;
    float relDev = 0.f;
    for (auto i = 0; i < roots.size(); ++i) {
        float minDevCur = numeric_limits<float>::max();
        auto tempInInd = 0;
        auto tempInd = 0;
        for (auto j = 0; j < roots.size(); ++j) {
            auto dev = abs(tempRootsInput[j] - tempRoots[i]);
            if (dev < minDevCur) {
                tempInInd = j;
                tempInd = i;
                minDevCur = dev;
            }
        }
        maxDev = max(minDevCur, maxDev);
        relDev = maxDev / abs(max(tempRootsInput[tempInInd], tempRoots[tempInd]));
        tempRootsInput[tempInInd] = -numeric_limits<float>::max();
        tempRoots[tempInd] = numeric_limits<float>::max();
    }
    return {maxDev, relDev};
}

/**
 * Generate root.
 *
 * @param low,high,dist values to generate root from low to high.
 * @return Root in based distance.
 */
void Framework::generate(float low, float high, float dist) {
    uniform_real_distribution<float> distribution(low, high);
    long double mid = distribution(generator);
    for (auto &root: roots) {
        uniform_real_distribution<float> distributionSmall(-dist, +dist);
        root = static_cast<float>(mid + distributionSmall(generator));
    }

    switch (roots.size()) {
        case 2: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1]);
            coefficients[2] = roots[0] * roots[1];
            break;
        }
        case 3: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1] + roots[2]);
            coefficients[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2];
            coefficients[3] = -(roots[0] * roots[1] * roots[2]);
            break;
        }
    }

    for (int i = 0; i < roots.size(); i++) {
        rootsReal[i] = static_cast<float>(roots[i].get_d());
    }
    for (int i = 0; i < roots.size() + 1; i++) {
        coefficientsReal[i] = static_cast<float>(coefficients[i].get_d());
    }
}


template<typename tmp>
void printVector(vector<tmp> input) {
    for (auto i = 0; i < input.size(); ++i) cout << input[i] << " ";
}

struct comparator {
    std::pair<float, float> deviations;
    vector<float> rootsCompute;
    vector<mpf_class> rootsTrue;
    vector<mpf_class> coefficients;
};

void Framework::generateBatch(int count, int rootsCount, float low, float high, float maxDistance,
                              vector<float> (*testing)(vector<float>)) {
    vector<Framework *> frameworks;

    float max_deviation = -1.f;
    auto cores = omp_get_num_procs();
    auto *comparison = new comparator[cores];

    for (auto i = 0; i < cores; ++i) {
        frameworks.push_back(new Framework(rootsCount));
        comparison[i].deviations = std::pair(-1.f, -1.f);
        comparison[i].rootsCompute = vector<float>(rootsCount);
        comparison[i].rootsTrue = vector<mpf_class>(rootsCount);
        comparison[i].coefficients = vector<mpf_class>(rootsCount + 1);
    }

#pragma OMP parallel for
    for (int i = 0; i < count; ++i) {
        //generation coefficients
        auto thread_id = omp_get_thread_num();
        frameworks[thread_id]->generate(low, high, maxDistance);

        //calculating roots
        auto outputRoots = testing(frameworks[thread_id]->coefficientsReal);
        auto deviation = frameworks[thread_id]->deviation(outputRoots);

        //calc deviation
        if (deviation.first > comparison[thread_id].deviations.first and
            deviation.first != numeric_limits<float>::infinity()) {
            comparison[thread_id].rootsCompute = vector<float>(outputRoots);
            comparison[thread_id].deviations = pair(deviation);
            comparison[thread_id].rootsTrue = vector<mpf_class>(frameworks[thread_id]->roots);
            comparison[thread_id].coefficients = vector<mpf_class>(frameworks[thread_id]->coefficients);
        }
    }

    //summ deviation
    comparator worse;
    for (auto i = 0; i < cores; ++i) {
        if (comparison[i].deviations.first > worse.deviations.first) {
            worse = comparison[i];
        }
    }

    cout << "Max deviation: " << worse.deviations.first << endl;
    cout << "Relative deviation: " << worse.deviations.second << endl;
    cout << "Bad computed roots: ";
    printVector(worse.rootsCompute);
    cout << endl;
    cout << "Original roots: ";
    printVector(worse.rootsTrue);
    cout << endl;
    cout << "Coefficients: ";
    printVector(worse.coefficients);
    cout << endl;

}

template float pr_product_difference<float>(float a, float b, float c, float d);

template double pr_product_difference<double>(double a, double b, double c, double d);

template long double pr_product_difference<long double>(long double a, long double b, long double c, long double d);

