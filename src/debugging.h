#pragma once
#include <iostream>

inline std::ostream &operator<<(std::ostream &out, const glm::ivec2 &ivec) {
    out << "(" << ivec.x << ", " << ivec.y << ")";
    return out;
}

inline std::ostream &operator<<(std::ostream &out, const glm::vec2 &vec) {
    out << "(" << vec.x << ", " << vec.y << ")";
    return out;
}

inline std::ostream &operator<<(std::ostream &out, const glm::vec3 &vec) {
    out << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return out;
}
