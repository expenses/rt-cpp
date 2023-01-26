#pragma once

struct Ray {
    Ray(RTCRayHit embree_ray);

    Ray(glm::vec3 o, glm::vec3 d, float max_t) : o(o), d(d), max_t(max_t) {
    }

    glm::vec3 o, d;
    float max_t;

    RTCRayHit as_embree();
};

struct Sphere {
    glm::vec3 center;
    float radius;
};
