#include "sampling.hpp"
#include <cstdio>
#include <random>
#include <span>
#include <vector>

int main() {
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<float> float_dist(0.0, 1.0);

    std::vector<float> values = {1.0f, 2.0f, 1.0f, 3.0f};
    std::vector<int> counts(values.size(), 0);
    /*auto dist = Distribution1D(values);


    for (int i = 0; i < 100000; i++) {
        auto index = dist.sample(float_dist(e1)).index;
        counts[index] ++;
    }

    for (int value: counts) {
        printf("%i, ",value);
    }
    printf("\n");*/

    std::vector<float> values_2d = {0.0, 1.0, 0.0, 1.0f};
    dbg(values_2d);

    int counts_2d[4] = {0, 0, 0, 0};

    Distribution2D dist_2d = Distribution2D(values_2d, 2, 2);

    /*
        for (int i = 0; i < 100000; i++) {
            auto indices = dist_2d.sample(vec2(float_dist(e1),
       float_dist(e1))).indices; counts_2d[indices.x + indices.y * 2] ++;
        }*/

    Sample2D samp = dist_2d.sample(vec2(0.75, 0.75));

    dbg(samp);
}
