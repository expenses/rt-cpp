#pragma once
#include "sampling.hpp"

struct Image {
    Image(const std::string &filename);

    const vec3 sample(glm::vec2 uv);

    const vec3 sample_env_map(glm::vec3 dir);

    const std::pair<vec3, vec3> sample_pdf(glm::vec2 rng);

    uint32_t width;
    uint32_t height;
    std::vector<vec3> rgb;
    Distribution2D distribution;
};
