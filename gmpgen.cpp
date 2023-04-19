#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "gmpgen.h"

using namespace std;


template<typename fp_t>
fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d) {
    auto tmp = d * c;
    return fma(a, b, -tmp) + fma(-d, c, tmp);
}

template<typename fp_t>
pair<fp_t, fp_t> Framework<fp_t>::deviation(vector<fp_t> rootsInput) {
    assert(5 > rootsInput.size() > 0);
    vector<fp_t> tempRoots = vector<fp_t>(roots.size());
    vector<fp_t> tempRootsInput = vector<fp_t>(roots.size());
    for (auto i = 0; i < roots.size(); ++i) {
        tempRoots[i] = roots[i];
        tempRootsInput[i] = rootsInput[i];
    }
    fp_t maxDev = 0;
    fp_t relDev = 0;
    for (auto i = 0; i < roots.size(); ++i) {
        fp_t minDevCur = numeric_limits<fp_t>::max();
        auto tempInInd = 0;
        auto tempInd = 0;
        for (auto j = 0; j < tempRootsInput.size(); ++j) {
            auto dev = abs(tempRootsInput[j] - tempRoots[i]);
            if (dev < minDevCur) {
                tempInInd = j;
                tempInd = i;
                minDevCur = dev;
            }
        }
        maxDev = max(minDevCur, maxDev);
        relDev = max(relDev, maxDev / max(abs(tempRootsInput[tempInInd]), abs(tempRoots[tempInd])));
        tempRootsInput[tempInInd] = -numeric_limits<fp_t>::max();
        tempRoots[tempInd] = numeric_limits<fp_t>::max();
    }
    return {maxDev, relDev};
}

template<typename fp_t>
void
Framework<fp_t>::generate(fp_t low, fp_t high, fp_t maxDistance, int multipleRoots, default_random_engine &generator) {
    uniform_real_distribution<fp_t> distribution(low, high);
    fp_t mid = distribution(generator);
    for (auto &root: roots) {
        uniform_real_distribution<fp_t> distributionSmall(-maxDistance, +maxDistance);
        root = mid + distributionSmall(generator);
    }
    for (auto i = 1; i < multipleRoots; ++i) {
        roots[i] = roots[0];
    }


    switch (roots.size()) {
        case 1: {
            coefficients[0] = 1;
            coefficients[1] = roots[0];
            break;
        }
        case 2: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1]);
            coefficients[2] = roots[0] * roots[1];
            break;
        }
        case 3: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1] + roots[2]);
//            coefficients[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2];
            coefficients[2] = fma(roots[1], roots[2], pr_product_difference(roots[0], roots[1], -roots[0], roots[2]));
            coefficients[3] = -(roots[0] * roots[1] * roots[2]);
            break;
        }
        case 4: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1] + roots[2] + roots[3]);
//            coefficients[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[0] * roots[3] + roots[1] * roots[2] +
//                              roots[1] * roots[3] + roots[2] * roots[3];
            coefficients[2] = pr_product_difference(roots[0], roots[1], -roots[0], roots[2]) +
                              pr_product_difference(roots[0], roots[3], -roots[1], roots[2]) +
                              pr_product_difference(roots[1], roots[3], -roots[2], roots[3]);
//            coefficients[3] = -(roots[0] * roots[1] * roots[2] + roots[0] * roots[1] * roots[3] +
//                                roots[0] * roots[2] * roots[3] + roots[1] * roots[2] * roots[3]);
            coefficients[3] = pr_product_difference(roots[0], fma(roots[2], roots[3],
                                                                  pr_product_difference(roots[1], roots[2], -roots[1],
                                                                                        roots[3])),
                                                    -roots[1] * roots[2], roots[3] * roots[3]);
            coefficients[4] = roots[0] * roots[1] * roots[2] * roots[3];
            break;
        }
    }
}

template<typename fp_t>
void Framework<fp_t>::generateSlow(fp_t low, fp_t high, fp_t maxDistance, int multipleRoots,
                                   default_random_engine &generator) {
    uniform_real_distribution<fp_t> distribution(low, high);
    fp_t mid = distribution(generator);
    vector<mpf_class> bigRoots(roots.size());
    vector<mpf_class> bigCoefficients(coefficients.size());
    for (auto &root: bigRoots) {
        uniform_real_distribution<fp_t> distributionSmall(-maxDistance, +maxDistance);
        root = mid + distributionSmall(generator);
    }
    for (auto i = 1; i < multipleRoots; ++i) {
        bigRoots[i] = bigRoots[0];
    }


    switch (bigRoots.size()) {
        case 1: {
            bigCoefficients[0] = 1;
            bigCoefficients[1] = bigRoots[0];
        }
        case 2: {
            bigCoefficients[0] = 1;
            bigCoefficients[1] = -(bigRoots[0] + bigRoots[1]);
            bigCoefficients[2] = bigRoots[0] * bigRoots[1];
            break;
        }
        case 3: {
            bigCoefficients[0] = 1;
            bigCoefficients[1] = -(bigRoots[0] + bigRoots[1] + bigRoots[2]);
            bigCoefficients[2] = bigRoots[0] * bigRoots[1] + bigRoots[0] * bigRoots[2] + bigRoots[1] * bigRoots[2];
            bigCoefficients[3] = -(bigRoots[0] * bigRoots[1] * bigRoots[2]);
            break;
        }
        case 4: {
            bigCoefficients[0] = 1;
            bigCoefficients[1] = -(bigRoots[0] + bigRoots[1] + bigRoots[2] + bigRoots[3]);
            bigCoefficients[2] = bigRoots[0] * bigRoots[1] + bigRoots[0] * bigRoots[2] + bigRoots[0] * bigRoots[3] +
                                 bigRoots[1] * bigRoots[2] +
                                 bigRoots[1] * bigRoots[3] + bigRoots[2] * bigRoots[3];
            bigCoefficients[3] = -(bigRoots[0] * bigRoots[1] * bigRoots[2] + bigRoots[0] * bigRoots[1] * bigRoots[3] +
                                   bigRoots[0] * bigRoots[2] * bigRoots[3] + bigRoots[1] * bigRoots[2] * bigRoots[3]);
            bigCoefficients[4] = bigRoots[0] * bigRoots[1] * bigRoots[2] * bigRoots[3];
            break;
        }
    }
    for (int i = 0; i < bigRoots.size(); ++i) {
        roots[i] = bigRoots[i].get_d();
    }
    for (int i = 0; i < bigCoefficients.size(); ++i) {
        coefficients[i] = bigCoefficients[i].get_d();
    }
}


template<typename tmp>
void printVector(vector<tmp> input) {
    for (auto i = 0; i < input.size(); ++i) cout << ' ' << input[i] << " " << endl;
}

template<typename fp_t>
struct comparator {
    std::pair<fp_t, fp_t> deviations;
    vector<fp_t> rootsCompute;
    vector<fp_t> rootsTrue;
    vector<fp_t> coefficients;
};

template<typename fp_t>
void Framework<fp_t>::generateBatch(int count, int rootsCount, fp_t low, fp_t high, fp_t maxDistance, int multipleRoots,
                                    vector<fp_t> (*testing)(vector<fp_t>), unsigned long long seed, bool slow) {
    assert(count > 0);
    assert(rootsCount > 0);
    assert(high - low > maxDistance > 0);
    assert(rootsCount >= multipleRoots >= 0);

    random_device randomDevice;
    auto cores = omp_get_num_procs();
    if (cores > count) cores = count;
    auto *generators = new default_random_engine[cores];
    auto *comparsion = new comparator<fp_t>[cores];
    auto *frameworks = new Framework<fp_t>[cores];

    if (!seed) seed = randomDevice();
#pragma OMP parallel for
    for (auto i = 0; i < cores; ++i) {
        frameworks[i] = Framework(rootsCount);
        generators[i] = default_random_engine(seed + i);
    }

    cout << "Started" << endl;
#pragma OMP parallel for
    for (auto i = 0; i < count; ++i) {
        //generation coefficients
        auto thread_id = omp_get_thread_num();
        if (slow) frameworks[thread_id].generateSlow(low, high, maxDistance, multipleRoots, generators[thread_id]);
        else frameworks[thread_id].generate(low, high, maxDistance, multipleRoots, generators[thread_id]);


        //calculating roots
        auto outputRoots = testing(frameworks[thread_id].coefficients);
        auto deviations = frameworks[thread_id].deviation(outputRoots);

        //calc deviation
        if (deviations.first >= comparsion[thread_id].deviations.first and
            deviations.first != numeric_limits<fp_t>::infinity()) {
            comparsion[thread_id].rootsCompute = vector<fp_t>(outputRoots);
            comparsion[thread_id].deviations = pair(deviations);
            comparsion[thread_id].rootsTrue = vector<fp_t>(frameworks[thread_id].roots);
            comparsion[thread_id].coefficients = vector<fp_t>(frameworks[thread_id].coefficients);
        }
    }

    //summ deviation
    comparator<fp_t> worse = comparsion[0];
    for (auto i = 0; i < cores; ++i) {
        if (comparsion[i].deviations.first > worse.deviations.first) {
            worse = comparsion[i];
        }
    }
    fixed(cout);
    cout.precision(numeric_limits<fp_t>::digits);
    cout << "Max deviation: " << worse.deviations.first << endl;
    cout << "Relative deviation: " << worse.deviations.second << endl;
    cout.precision(numeric_limits<fp_t>::digits10);
    cout << "Bad computed roots: " << endl;
    printVector(worse.rootsCompute);
    cout << "Original roots: " << endl;
    printVector(worse.rootsTrue);
    cout << "Coefficients: " << endl;
    printVector(worse.coefficients);

    delete[] frameworks;
    delete[] comparsion;
    delete[] generators;
}

template float pr_product_difference<float>(float a, float b, float c, float d);

template double pr_product_difference<double>(double a, double b, double c, double d);

template long double pr_product_difference<long double>(long double a, long double b, long double c, long double d);

template void
Framework<float>::generateBatch(int count, int rootsCount, float low, float high, float maxDistance, int cRoots,
                                vector<float> (*testing)(vector<float>), unsigned long long seed = 0,
                                bool slow = false);

template void
Framework<double>::generateBatch(int count, int rootsCount, double low, double high, double maxDistance, int cRoots,
                                 vector<double> (*testing)(vector<double>), unsigned long long seed = 0,
                                 bool slow = false);

