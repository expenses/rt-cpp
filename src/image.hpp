#pragma once
#include "sampling.hpp"

struct Image {
    Image(const std::string &filename, OIIO::TextureSystem &texture_system);

    const vec3 sample(glm::vec2 uv);

    const vec3 sample_env_map(glm::vec3 dir);

    const std::tuple<float, vec2, vec3> sample_pdf(glm::vec2 rng);

    Distribution2D distribution;
    OIIO::ustring path;
    OIIO::TextureSystem &texture_system;
};
