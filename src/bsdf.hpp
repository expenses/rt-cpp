#pragma once
#include "external.hpp"
#include <cmath>

using namespace glm;

inline float cos_theta(vec3 w) {
    return w.z;
}

inline float cos_2_theta(vec3 w) {
    return cos_theta(w) * cos_theta(w);
}

inline float sin_2_theta(vec3 w) {
    return max(0.0f, 1.0f - cos_2_theta(w));
}

inline float sin_theta(vec3 w) {
    return sqrt(sin_2_theta(w));
}

inline float tan_theta(vec3 w) {
    return sin_theta(w) / cos_theta(w);
}

inline float tan_2_theta(vec3 w) {
    return sin_2_theta(w) / cos_2_theta(w);
}

inline float cos_phi(vec3 w) {
    auto sin_theta_ = sin_theta(w);
    return (sin_theta_ == 0) ? 1 : clamp(w.x / sin_theta_, -1.0f, 1.0f);
}

inline float sin_phi(vec3 w) {
    auto sin_theta_ = sin_theta(w);
    return (sin_theta_ == 0) ? 0 : clamp(w.y / sin_theta_, -1.0f, 1.0f);
}

inline float cos_2_phi(vec3 w) {
    return cos_phi(w) * cos_phi(w);
}

inline float sin_2_phi(vec3 w) {
    return sin_phi(w) * sin_phi(w);
}

struct TrowbridgeReitzDistribution {
    vec2 alpha;

    float eval(vec3 normal) {
        float tan_2_theta_ = tan_2_theta(normal);

        if (std::isinf(tan_2_theta_)) {
            return 0;
        }

        auto cos_4_theta = cos_2_theta(normal) * cos_2_theta(normal);

        auto alpha_2 = alpha * alpha;

        auto e = (cos_2_phi(normal) / alpha_2.x + sin_2_phi(normal) / alpha_2.y) * tan_2_theta_;

        return 1.0f / (M_PI * alpha.x * alpha.y * cos_4_theta * (1 + e) * (1 + e));
    }
};

struct MicrofacetReflection {
    TrowbridgeReitzDistribution distribution;
};

// Use a tagged enum beacuse std::variant uses vtables :(
struct Bsdf {
    enum Tag
    {
        Diffuse,
        Conductor,
    } tag;

    union Params {
        struct Diffuse {
            vec3 colour;
        } diffuse;
        struct Conductor {
            vec3 colour;
            vec2 alpha;
        } conductor;
    } params;
};
