#include <chrono>
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <iostream>
#include <math.h>
#include <random>

#include "buffers.cpp"
#include "fexr.cpp"
#include "geo.cpp"
#include "sampling.hpp"
#include "tev.cpp"
#include "util.hpp"

using namespace glm;

int main() {
    std::string ff = "fexr/san.fexr"; //"neon_photostudio_4k.exr.fexr";
    auto env_map = Fexr(ff);

    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<float> float_dist(0.0, 1.0);

    TevConnection connection = TevConnection();

    printf("Started\n");

    std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    uint32_t width = 512;
    uint32_t height = 512;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    auto start = std::chrono::steady_clock::now();
    bool never_updated = true;

    auto red = vec3(1.0, 0.1, 0.1);
    auto green = vec3(0.2, 0.9, 0.2);
    auto blue = vec3(0.3, 0.3, 0.8);

    Sphere spheres[8] = {
        Sphere{vec3(-0.25, -0.25, -0.25), 0.25, red},
        Sphere{vec3(-0.25, +0.25, -0.25), 0.15, green},
        Sphere{vec3(-0.25, -0.25, +0.25), 0.15, blue},
        Sphere{vec3(-0.25, +0.25, +0.25), 0.15, red},
        Sphere{vec3(+0.25, -0.25, -0.25), 0.25, green},
        Sphere{vec3(+0.25, +0.25, -0.25), 0.25, blue},
        Sphere{vec3(+0.25, -0.25, +0.25), 0.25, green},
        Sphere{vec3(+0.25, +0.25, +0.25), 0.15, blue},

        // Sphere{vec3(0, -100.5, 0), 100, vec3(0.1, 1.0, 0.1)},
    };

    vec3 origin = vec3(-1.0, 0.5, 0.0);
    vec3 look_at = vec3(0.0, 0.0, 0.0);

    mat4 view = lookAt(origin, look_at, vec3(0.0, 1.0, 0.0));
    mat4 proj = perspective(radians(75.0f), (float)width / (float)height,
                            0.001f, 10000.0f);

    mat4 view_inverse = inverse(view);
    mat4 proj_inverse = inverse(proj);

    while (true) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float u = ((float)x + float_dist(e1)) / width;
                float v = ((float)y + float_dist(e1)) / height;

                auto target = proj_inverse *
                              vec4(u * 2.0 - 1.0, v * -2.0 + 1.0, 1.0, 1.0);
                auto local_direction = normalize(vec3(target));
                auto direction =
                    normalize(vec3(view_inverse * vec4(local_direction, 0.0)));

                auto ray = Ray{origin, direction, 10000.0f};

                Intersection closest_intersection;
                closest_intersection.valid = false;
                vec3 sphere_colour;

                for (auto &sphere : spheres) {
                    auto intersection = sphere.intersect(ray);

                    if (intersection.valid &&
                        (!closest_intersection.valid ||
                         intersection.t < closest_intersection.t)) {
                        closest_intersection = intersection;
                        sphere_colour = sphere.colour;
                        ray.max_t = intersection.t;
                    }
                }

                vec3 colour;

                if (closest_intersection.valid) {
                    auto pdf_sample = env_map.sample_pdf(
                        vec2(float_dist(e1), float_dist(e1)));

                    colour = pdf_sample.value * sphere_colour *
                             vec3(std::max(dot(pdf_sample.direction,
                                               closest_intersection.normal),
                                           0.0f)) /
                             vec3(M_PI); // env_map.sample_env_map(dir) *
                                         // sphere.colour;

                    auto shadow_ray =
                        Ray{closest_intersection.position +
                                pdf_sample.direction * vec3(0.001),
                            pdf_sample.direction, 10000.0f};

                    for (auto &sphere : spheres) {
                        if (sphere.intersect(shadow_ray).valid) {
                            colour = vec3(0.0);
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
            printf("accum.num_samples: %lu\n", accum.num_samples);

            start = end;
            never_updated = false;

            accum.update_output(output);

            connection.send_update(image_name, 0, 0, width, height,
                                   channel_names, channel_offsets,
                                   channel_strides, output.data);
        }
    }
}
