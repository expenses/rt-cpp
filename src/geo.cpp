#include <glm/geometric.hpp>
#include <glm/vec3.hpp>
#include <utility>

#include "util.hpp"

using namespace glm;

struct Intersection {
    float t;
    vec3 normal;
    vec3 position;
    bool valid;
};

struct Ray {
    vec3 o, d;
    float max_t;
};

struct Sphere {
    vec3 center;
    float radius;
    vec3 colour;

    Intersection intersect(const Ray &ray) {
        auto oc = ray.o - center;

        auto a = dot(ray.d, ray.d);
        auto b = 2.0 * dot(oc, ray.d);
        auto c = dot(oc, oc) - radius * radius;

        Intersection intersection = {0.0, vec3(0.0), vec3(0.0), false};

        auto roots = quadratic(a, b, c);

        if (!roots.valid) {
            return intersection;
        }

        if (roots.t_1 <= 0.0 || roots.t_0 > ray.max_t) {
            return intersection;
        }

        intersection.valid = true;
        intersection.t = roots.t_0;

        if (roots.t_0 <= 0.0) {
            intersection.t = roots.t_1;
        }

        intersection.position = ray.o + intersection.t * ray.d;
        intersection.normal = normalize(intersection.position - center);

        return intersection;
    }
};
