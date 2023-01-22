#pragma once
#include "geo.h"
#include "image.h"

struct Scene {
    std::vector<Sphere> spheres;
    Image environment_map;

    std::optional<std::tuple<Intersection, Bsdf>> find_intersection(Ray ray, bool find_any);
};
