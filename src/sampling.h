#pragma once

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <span>
#include <vector>

using namespace glm;

struct Sample1D {
    const int index;
    const float pdf;
    const float point;
};

std::ostream &operator<<(std::ostream &out, const Sample1D &sample);

struct Distribution1D {
    Distribution1D();

    Distribution1D(std::span<float> func);

    Sample1D sample(float rng);

    std::vector<float> cdf;
    std::vector<float> function;
    float function_integral;
};

struct Sample2D {
    float pdf;
    ivec2 indices;
    vec2 point;
};

std::ostream &operator<<(std::ostream &out, const Sample2D &sample);

struct Distribution2D {
    Distribution2D();

    Distribution2D(std::span<float> image, uint32_t width, uint32_t height);

    Sample2D sample(vec2 rng);

    std::vector<Distribution1D> row_distributions;
    Distribution1D main_distribution;
};
