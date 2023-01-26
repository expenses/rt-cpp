#pragma once
#include "geo.hpp"
#include "image.hpp"
#include "bsdf.hpp"

struct Scene {
    std::vector<Bsdf> sphere_bsdfs;
    Image environment_map;
};
