#pragma once
#include "bsdf.hpp"
#include "geo.hpp"
#include "image.hpp"

struct Scene {
    std::vector<Bsdf> sphere_bsdfs;
    Image environment_map;
};
