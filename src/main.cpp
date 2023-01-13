#include <chrono>
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <iostream>
#include <cmath>
#include <random>

#include "buffers.cpp"
#include "geo.cpp"
#include "image.cpp"
#include "sampling.hpp"
#include "tev.cpp"
#include "util.hpp"

using namespace glm;

std::optional<std::tuple<Intersection, Bsdf>> find_intersection(std::span<Sphere> spheres, Ray ray, bool find_any) {
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

int main() {
    auto env_map = Image("san.exr");

    std::random_device rand_dev;
    std::default_random_engine rng(rand_dev());
    std::uniform_real_distribution<float> float_dist(0.0, 1.0);

    TevConnection connection = TevConnection();

    const std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    const uint32_t width = 512;
    const uint32_t height = width;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    auto green = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.2, 0.9, 0.2)}}};
    auto blue = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.3, 0.3, 0.8)}}};

    auto diffuse = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.5, 0.5, 0.25)}}};

    auto conductor_a =
        Bsdf{.tag = Bsdf::Tag::Conductor,
             .params = Bsdf::Params{.conductor = {.colour = vec3(0.5, 0.5, 0.25), .alpha = vec2(0.9, 0.001)}}};

    auto conductor_b =
        Bsdf{.tag = Bsdf::Tag::Conductor,
             .params = Bsdf::Params{.conductor = {.colour = vec3(0.5, 0.5, 0.25), .alpha = vec2(0.001, 0.9)}}};

    auto conductor_c =
        Bsdf{.tag = Bsdf::Tag::Conductor,
             .params = Bsdf::Params{.conductor = {.colour = vec3(0.5, 0.5, 0.25), .alpha = vec2(0.15, 0.15)}}};

    auto conductor_d =
        Bsdf{.tag = Bsdf::Tag::Conductor,
             .params = Bsdf::Params{.conductor = {.colour = vec3(0.5, 0.5, 0.25), .alpha = vec2(0.03, 0.03)}}};

    Sphere spheres[4] = {
        // Sphere{.center = vec3(-0.25, -0.25, -0.25), .radius = 0.25, .bsdf = red},
        // Sphere{.center = vec3(-0.25, +0.25, -0.25), .radius = 0.15, .bsdf = green},
        // Sphere{.center = vec3(-0.25, -0.25, +0.25), .radius = 0.15, .bsdf = blue},
        // Sphere{.center = vec3(-0.25, +0.25, +0.25), .radius = 0.15, .bsdf = red},
        // Sphere{.center = vec3(+0.25, -0.25, -0.25), .radius = 0.25, .bsdf = green},
        // Sphere{.center = vec3(+0.25, +0.25, -0.25), .radius = 0.25, .bsdf = blue},
        Sphere{.center = vec3(0, -0.6, 0), .radius = 0.5, .bsdf = conductor_a},
        Sphere{.center = vec3(0, 0.5, 0), .radius = 0.5, .bsdf = conductor_b},
        Sphere{.center = vec3(0.75, 0.5, 0.75), .radius = 0.5, .bsdf = conductor_d},
        Sphere{.center = vec3(0.75, -0.6, 0.75), .radius = 0.5, .bsdf = diffuse},
        // Sphere{.center = vec3(+0.25, +0.25, +0.25), .radius = 0.15, .bsdf = blue},

        // Sphere{vec3(0, -100.5, 0), 100, vec3(0.1, 1.0, 0.1)},
    };

    const vec3 origin = vec3(1.5, 0.2, -1.5);
    const vec3 look_at = vec3(0.0, 0.0, 0.75);

    const mat4 view = lookAt(origin, look_at, vec3(0.0, 1.0, 0.0));
    const mat4 proj = perspective(radians(59.0f), float(width) / float(height), 0.001f, 10000.0f);

    const mat4 view_inverse = inverse(view);
    const mat4 proj_inverse = inverse(proj);

    auto start = std::chrono::steady_clock::now();
    bool never_updated = true;

    while (true) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                const float u = (float(x) + float_dist(rng)) / width;
                const float v = (float(y) + float_dist(rng)) / height;

                auto target = proj_inverse * vec4(u * 2.0 - 1.0, v * -2.0 + 1.0, 1.0, 1.0);
                auto local_direction = normalize(vec3(target));
                auto direction = normalize(vec3(view_inverse * vec4(local_direction, 0.0)));

                auto ray = Ray{origin, direction, 10000.0f};

                vec3 colour = vec3(0.0);

                if (auto intersection_ = find_intersection(spheres, ray, false)) {
                    auto [intersection, bsdf] = intersection_.value();

                    switch (bsdf.tag) {
                    case Bsdf::Tag::Diffuse: {
                        auto params = bsdf.params.diffuse;

                        auto [value, direction] = env_map.sample_pdf(vec2(float_dist(rng), float_dist(rng)));

                        auto shadow_ray = Ray{intersection.position + direction * vec3(0.0001), direction, 10000.0f};

                        if (!find_intersection(spheres, shadow_ray, true)) {
                            auto l_dot_n = std::max(dot(direction, intersection.normal), 0.0f);
                            colour = value * params.colour * vec3(l_dot_n / M_PI);
                        }

                        break;
                    }
                    case Bsdf::Tag::Conductor: {
                        auto params = bsdf.params.conductor;

                        auto hemisphere_sample = CosineSampleHemisphere(vec2(float_dist(rng), float_dist(rng)));
                        hemisphere_sample.x *= params.alpha.x;
                        hemisphere_sample.z *= params.alpha.y;
                        hemisphere_sample = normalize(hemisphere_sample);

                        auto sample_dir = rotation_matrix(reflect(ray.d, intersection.normal)) * hemisphere_sample;

                        auto reflected_ray =
                            Ray{intersection.position + sample_dir * vec3(0.0001), sample_dir, 10000.0f};

                        if (!find_intersection(spheres, reflected_ray, true)) {
                            colour = env_map.sample_env_map(reflected_ray.d) * params.colour;
                        }

                        break;
                    }
                    }
                } else {
                    colour = env_map.sample_env_map(ray.d);
                }

                auto offset = (y * width + x) * 3;
                accum.data[offset + 0] += colour.x;
                accum.data[offset + 1] += colour.y;
                accum.data[offset + 2] += colour.z;
            }
        }

        accum.num_samples += 1;

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        if (elapsed_seconds.count() > 1.0 || never_updated) {
            dbg(accum.num_samples);

            start = end;
            never_updated = false;

            accum.update_output(output);

            connection.send_update(image_name, 0, 0, width, height, channel_names, channel_offsets, channel_strides,
                                   output.data);
        }
    }
}
