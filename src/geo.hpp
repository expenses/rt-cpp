#pragma once

using namespace glm;

struct Sphere {
    vec3 center;
    float radius;
};

struct alignas(16) Ray {
    vec3 origin;
    float t_near;
    vec3 direction;
    float time;
    float t_far;
    uint mask;
    uint id;
    uint flags;

    vec3 at(float t) {
        return origin + direction * t;
    }
};

struct alignas(16) Hit {
    vec3 normal;
    vec2 uv;
    uint primID;
    uint geomID;
    uint instID[RTC_MAX_INSTANCE_LEVEL_COUNT];
};

struct alignas(16) RayHit {
    Ray ray;
    Hit hit;
};
