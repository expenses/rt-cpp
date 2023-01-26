#include "geo.hpp"

#include <glm/geometric.hpp>
#include <glm/vec3.hpp>
#include <optional>
#include <utility>

#include "util.hpp"

Ray::Ray(RTCRayHit embree_ray) {
    o = vec3(embree_ray.ray.org_x, embree_ray.ray.org_y, embree_ray.ray.org_z);
    d = vec3(embree_ray.ray.dir_x, embree_ray.ray.dir_y, embree_ray.ray.dir_z);
    max_t = embree_ray.ray.tfar;
}

RTCRayHit Ray::as_embree() {
    return RTCRayHit{.ray =
                         {
                             .org_x = o.x,
                             .org_y = o.y,
                             .org_z = o.z,
                             .tnear = 0.0,
                             .dir_x = d.x,
                             .dir_y = d.y,
                             .dir_z = d.z,
                             .tfar = max_t,
                         },
                     .hit = {.geomID = RTC_INVALID_GEOMETRY_ID}};
}
