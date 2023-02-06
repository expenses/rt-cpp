#pragma once

#include <glm/geometric.hpp>
#include <glm/mat3x3.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

using namespace glm;


inline vec3 pxr_to_vec3(pxr::GfVec3f vec) {
    return vec3(vec[0], vec[1], vec[2]);
}

inline pxr::GfVec3f vec3_to_pxr(vec3 vec) {
    return pxr::GfVec3f(vec.x, vec.y, vec.z);
}


inline vec2 UniformSampleDisk(vec2 u) {
    auto r = sqrt(u.x);
    auto theta = 2.0 * M_PI * u.y;
    return vec2(r * cos(theta), r * sin(theta));
}

inline vec3 CosineSampleHemisphere(vec2 rng) {
    auto d = UniformSampleDisk(rng);
    auto z = sqrt(1.0 - d.x * d.x - d.y * d.y);
    if (std::isnan(z)) {
        z = 0.0;
    }
    return vec3(d.x, z, d.y);
}

inline mat3 rotation_matrix(vec3 new_y) {
    if (new_y.y == 1.0) {
        return mat3(vec3(1,0,0), new_y, vec3(0,0,1));
    }

    vec3 new_z = normalize(cross(new_y, vec3(0, 1, 0)));
    vec3 new_x = normalize(cross(new_y, new_z));
    return mat3(new_x, new_y, new_z);
}

inline float luminance(vec3 rgb) {
    return rgb.r * 0.212671f + rgb.g * 0.715160f + rgb.b * 0.072169f;
}
