#pragma once
#include "bsdf.hpp"
#include <optional>

struct Intersection {
    float t;
    glm::vec3 normal;
    glm::vec3 position;
};

struct Ray {
    glm::vec3 o, d;
    float max_t;
};

struct Sphere {
    glm::vec3 center;
    float radius;
    Bsdf bsdf;

    std::optional<Intersection> intersect(const Ray &ray);
};
