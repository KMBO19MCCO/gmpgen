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
        tempRoots[i] = static_cast<fp_t>(roots[i].get_d());
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
void Framework<fp_t>::generate(fp_t low, fp_t high, fp_t maxDistance, int multipleRoots, mt19937_64 generator) {


    uniform_real_distribution<fp_t> distribution(low, high);
    long double mid = distribution(generator);
    for (auto &root: roots) {
        uniform_real_distribution<fp_t> distributionSmall(-maxDistance, +maxDistance);
        root = static_cast<fp_t>(mid + distributionSmall(generator));
    }
    for (auto i = 1; i < multipleRoots; ++i) {
        roots[i] = roots[0];
    }

    switch (roots.size()) {
        case 1: {
            coefficients[0] = 1;
            coefficients[1] = roots[0];
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
            coefficients[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2];
            coefficients[3] = -(roots[0] * roots[1] * roots[2]);
            break;
        }
        case 4: {
            coefficients[0] = 1;
            coefficients[1] = -(roots[0] + roots[1] + roots[2] + roots[3]);
            coefficients[2] = roots[0] * roots[1] + roots[0] * roots[2] + roots[0] * roots[3] + roots[1] * roots[2] +
                              roots[1] * roots[3] + roots[2] * roots[3];
            coefficients[3] = -(roots[0] * roots[1] * roots[2] + roots[0] * roots[1] * roots[3] +
                                roots[0] * roots[2] * roots[3] + roots[1] * roots[2] * roots[3]);
            coefficients[4] = roots[0] * roots[1] * roots[2] * roots[3];
            break;
        }
    }

    for (auto i = 0; i < roots.size(); ++i) {
        rootsReal[i] = static_cast<fp_t>(roots[i].get_d());
    }
    for (auto i = 0; i < roots.size() + 1; ++i) {
        coefficientsReal[i] = static_cast<fp_t>(coefficients[i].get_d());
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
    vector<mpf_class> rootsTrue;
    vector<mpf_class> coefficients;
};

template<typename fp_t>
void Framework<fp_t>::generateBatch(int count, int rootsCount, fp_t low, fp_t high, fp_t maxDistance, int multipleRoots,
                                    vector<fp_t> (*testing)(vector<fp_t>), unsigned long long seed) {
    assert(count > 0);
    assert(rootsCount > 0);
    assert(high - low > maxDistance > 0);
    assert(rootsCount >= multipleRoots >= 0);
    vector<Framework *> frameworks;

    auto cores = omp_get_num_procs();
    auto *comparison = new comparator<fp_t>[cores];

    for (auto i = 0; i < cores; ++i) {
        frameworks.push_back(new Framework(rootsCount));
        comparison[i].deviations = std::pair(static_cast<fp_t>(-1), static_cast<fp_t>(-1));
        comparison[i].rootsCompute = vector<fp_t>(rootsCount);
        comparison[i].rootsTrue = vector<mpf_class>(rootsCount);
        comparison[i].coefficients = vector<mpf_class>(rootsCount + 1);
    }

    random_device randomDevice;
    auto *generators = new mt19937_64[cores];
    if (!seed) seed = randomDevice();
    for (auto i = 0; i < cores; ++i) {
        generators[i] =  mt19937_64(seed);
    }

#pragma OMP parallel for
    for (auto i = 0; i < count; ++i) {
        //generation coefficients
        auto thread_id = omp_get_thread_num();
        frameworks[thread_id]->generate(low, high, maxDistance, multipleRoots, generators[thread_id]);

        //calculating roots
        auto outputRoots = testing(frameworks[thread_id]->coefficientsReal);
        auto deviations = frameworks[thread_id]->deviation(outputRoots);

        //calc deviation
        if (deviations.first > comparison[thread_id].deviations.first and
            deviations.first != numeric_limits<fp_t>::infinity()) {
            comparison[thread_id].rootsCompute = vector<fp_t>(outputRoots);
            comparison[thread_id].deviations = pair(deviations);
            comparison[thread_id].rootsTrue = vector<mpf_class>(frameworks[thread_id]->roots);
            comparison[thread_id].coefficients = vector<mpf_class>(frameworks[thread_id]->coefficients);
        }
    }

    //summ deviation
    comparator<fp_t> worse;
    for (auto i = 0; i < cores; ++i) {
        if (comparison[i].deviations.first > worse.deviations.first) {
            worse = comparison[i];
        }
    }

    cout << "Max deviation: " << worse.deviations.first << endl;
    cout << "Relative deviation: " << worse.deviations.second << endl;
    cout << "Bad computed roots: " << endl;
    printVector(worse.rootsCompute);
    cout << "Original roots: " << endl;
    printVector(worse.rootsTrue);
    cout << "Coefficients: " << endl;
    printVector(worse.coefficients);
}

template float pr_product_difference<float>(float a, float b, float c, float d);

template double pr_product_difference<double>(double a, double b, double c, double d);

template long double pr_product_difference<long double>(long double a, long double b, long double c, long double d);

template void
Framework<float>::generateBatch(int count, int rootsCount, float low, float high, float maxDistance, int cRoots,
                                vector<float> (*testing)(vector<float>),unsigned long long seed = 0);

template void
Framework<double>::generateBatch(int count, int rootsCount, double low, double high, double maxDistance, int cRoots,
                                 vector<double> (*testing)(vector<double>), unsigned long long seed = 0);

