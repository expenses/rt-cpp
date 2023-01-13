#include <glm/geometric.hpp>
#include <glm/vec3.hpp>
#include <optional>
#include <utility>

#include "bsdf.hpp"
#include "util.hpp"

using namespace glm;

struct Intersection {
    float t;
    vec3 normal;
    vec3 position;
};

struct Ray {
    vec3 o, d;
    float max_t;
};

struct Sphere {
    vec3 center;
    float radius;
    Bsdf bsdf;

    std::optional<Intersection> intersect(const Ray &ray) {
        auto oc = ray.o - center;

        auto a = dot(ray.d, ray.d);
        auto b = 2.0 * dot(oc, ray.d);
        auto c = dot(oc, oc) - radius * radius;

        auto roots = quadratic(a, b, c);

        if (!roots) {
            return std::nullopt;
        }

        auto [t_0, t_1] = roots.value();

        if (t_1 <= 0.0 || t_0 > ray.max_t) {
            return std::nullopt;
        }

        auto t = t_0;

        if (t_0 <= 0.0) {
            t = t_1;
        }

        auto position = ray.o + t * ray.d;

        return Intersection{
            .t = t,
            .normal = normalize(position - center),
            .position = position,
        };
    }
};
