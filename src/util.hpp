#pragma once

#include <glm/geometric.hpp>
#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <optional>
#include <tuple>

#include "util.hpp"

using namespace glm;

vec2 UniformSampleDisk(vec2 u) {
    auto r = sqrt(u.x);
    auto theta = 2.0 * M_PI * u.y;
    return vec2(r * cos(theta), r * sin(theta));
}

vec3 CosineSampleHemisphere(vec2 rng) {
    auto d = UniformSampleDisk(rng);
    auto z = sqrt(1.0 - d.x * d.x - d.y * d.y);
    return vec3(d.x, z, d.y);
}

mat3 rotation_matrix(vec3 new_y) {
    vec3 new_z = normalize(cross(new_y, vec3(0, 1, 0)));
    vec3 new_x = normalize(cross(new_y, new_z));
    return mat3(new_x, new_y, new_z);
}

std::optional<std::tuple<float, float>> quadratic(double a, double b, double c) {
    auto discrim = b * b - 4.0 * a * c;

    if (discrim < 0.0) {
        return std::nullopt;
    }

    auto root_discrim = sqrt(discrim);

    double q = -0.5 * (b + root_discrim);

    if (b < 0.0) {
        q = -0.5 * (b - root_discrim);
    }

    auto t_0 = static_cast<float>(q / a);
    auto t_1 = static_cast<float>(c / q);

    if (t_0 > t_1) {
        std::swap(t_0, t_1);
    }

    return std::make_tuple(t_0, t_1);
}

float luminance(vec3 rgb) {
    return rgb.r * 0.212671f + rgb.g * 0.715160f + rgb.b * 0.072169f;
}
