#include <iostream>
#include <gmpxx.h>
#include <random>

using namespace std;

#define MAX_DISTANCE 1e-4

int main() {
    mpf_class a, b, c;
    mpf_class x0, x1;

    random_device randomDevice;
    mt19937_64 generator(randomDevice());
//    mt19937_64 generator(1337);
    uniform_real_distribution<float> distribution(-1.0, 1.0);
    long double mid = distribution(generator);
    uniform_real_distribution<float> distributionSmall(-MAX_DISTANCE, +MAX_DISTANCE);

    x0 = static_cast<float>(mid + distributionSmall(generator));
    x1 = static_cast<float>(mid + distributionSmall(generator));

    x0 = "0.7000558376312255859375";
    x1 = "0.6999988555908203125";

    cout.precision(50);
    cout << x0 << " " << x1 << endl;

    a = 1;
    b = -(x0 + x1);
    c = x0 * x1;
    cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
    return 0;
}