#include "sampling.hpp"
#include "debugging.hpp"
#include <cstdint>
#include <cstdio>

#include <span>
#include <vector>

using namespace glm;

std::ostream &operator<<(std::ostream &out, const Sample1D &sample) {
    out << "{pdf: " << sample.pdf << ", indices: " << sample.index << ", point: " << sample.point << "}";
    return out;
}

Distribution1D::Distribution1D(std::span<float> func) {
    function.reserve(func.size());
    function.insert(function.begin(), func.begin(), func.end());

    const auto n = function.size();

    cdf.resize(function.size() + 1);
    cdf[0] = 0.0;

    // The function's domain is split over n pieces so we divide by that:
    // https://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Sampling_Random_Variables#x1-Example:Piecewise-Constant1DFunctions.
    // This also results in the integral being the average of each piece.
    for (int i = 1; i < n + 1; i++) {
        cdf[i] = cdf[i - 1] + function[i - 1] / float(n);
    }

    function_integral = cdf[function.size()];

    for (float &value : cdf) {
        value /= function_integral;
    }
}

Sample1D Distribution1D::sample(float rng) {
    const auto pointer = lower_bound(cdf.begin(), cdf.end(), rng);
    const int upper_index = std::max(std::distance(cdf.begin(), pointer), 1l);
    const int lower_index = upper_index - 1;

    // Because the integral is equivalent of the average of each piece of
    // the function, the pdf will be >1 for higher-than-average values and
    // <1 for lower-than-average values.
    const auto pdf = function[lower_index] / function_integral;

    // Calculate the offset between the two indices as a fraction
    float offset_from_index = rng - cdf[lower_index];
    offset_from_index /= (cdf[upper_index] - cdf[lower_index]);

    const auto point = (float(lower_index) + offset_from_index) / function.size();

    return Sample1D{.index = lower_index, .pdf = pdf, .point = point};
}

std::ostream &operator<<(std::ostream &out, const Sample2D &sample) {
    out << "{pdf: " << sample.pdf << ", indices: " << sample.indices << ", point: " << sample.point << "}";
    return out;
}

Distribution2D::Distribution2D(std::span<float> image, uint32_t width, uint32_t height) {
    row_distributions.reserve(height);
    std::vector<float> row_integrals(height, 0.0);

    for (int y = 0; y < height; y++) {
        std::span<float> row_span = {image.data() + y * width, width};
        row_distributions.push_back(Distribution1D(row_span));
        row_integrals[y] = row_distributions[y].function_integral;
    }

    main_distribution = Distribution1D(row_integrals);
};

Sample2D Distribution2D::sample(vec2 rng) {
    Sample1D row_sample = main_distribution.sample(rng.x);
    Sample1D col_sample = row_distributions[row_sample.index].sample(rng.y);

    return Sample2D{.pdf = row_sample.pdf * col_sample.pdf,
                    .indices = ivec2(col_sample.index, row_sample.index),
                    .point = vec2(col_sample.point, row_sample.point)};
}
