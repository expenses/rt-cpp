#pragma once
#include <glm/glm.hpp>
#include <iostream>

using namespace glm;

std::ostream &operator<<(std::ostream &out, const ivec2 &ivec) {
    out << "(" << ivec.x << ", " << ivec.y << ")";
    return out;
}

std::ostream &operator<<(std::ostream &out, const vec2 &vec) {
    out << "(" << vec.x << ", " << vec.y << ")";
    return out;
}

std::ostream &operator<<(std::ostream &out, const vec3 &vec) {
    out << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return out;
}
