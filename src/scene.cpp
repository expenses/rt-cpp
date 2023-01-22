#include "scene.h"

std::optional<std::tuple<Intersection, Bsdf>> Scene::find_intersection(Ray ray, bool find_any) {
    std::optional<std::tuple<Intersection, Bsdf>> closest = std::nullopt;

    for (auto &sphere : spheres) {
        if (auto intersection_ = sphere.intersect(ray)) {
            auto intersection = intersection_.value();

            if (!closest || intersection.t < std::get<0>(closest.value()).t) {
                closest = std::optional(std::tuple(intersection, sphere.bsdf));
                ray.max_t = intersection.t;
            }

            if (find_any) {
                return closest;
            }
        }
    }

    return closest;
}
