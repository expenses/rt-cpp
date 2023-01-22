#include <cmath>
#include <iostream>

#include "external.hpp"

#include "buffers.hpp"
#include "debugging.hpp"
#include "geo.hpp"
#include "image.hpp"
#include "sampling.hpp"
#include "scene.hpp"
#include "tev.hpp"
#include "util.hpp"
#include <pxr/pxr.h>

using namespace glm;

int main(int argc, char *argv[]) {
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(argv[1]);
    pxr::UsdGeomXformCache cache;

    std::vector<Sphere> spheres;
    vec3 origin = vec3(0.0, 0.0, 0.0);
    vec3 look_at = vec3(0.0, 0.0, 0.0);
    float fov = 0.0;

    for (pxr::UsdPrim prim : pxr::UsdPrimRange::Stage(stage)) {
        if (prim.IsA<pxr::UsdGeomSphere>()) {
            auto sphere = pxr::UsdGeomSphere(prim);

            double radius;
            sphere.GetRadiusAttr().Get(&radius);

            auto transform = cache.GetLocalToWorldTransform(prim);
            pxr::GfVec3d translation = transform.ExtractTranslation();

            pxr::VtArray<pxr::GfVec3f> display_colours;
            prim.GetAttribute(pxr::UsdGeomTokensType().primvarsDisplayColor).Get(&display_colours);
            vec3 col = vec3(0.5, 0.5, 0.5);
            if (display_colours.size() > 0) {
                col = vec3(display_colours[0][0], display_colours[0][1], display_colours[0][2]);
            }

            auto diffuse = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = col}}};

            spheres.push_back(Sphere{.center = vec3(translation[0], translation[1], translation[2]),
                                     .radius = float(radius),
                                     .bsdf = diffuse});
        } else if (prim.IsA<pxr::UsdGeomCamera>()) {
            auto camera = pxr::UsdGeomCamera(prim).GetCamera(pxr::UsdTimeCode());
            fov = camera.GetFieldOfView(pxr::GfCamera::FOVDirection::FOVHorizontal);
            auto transform = camera.GetTransform();
            pxr::GfVec3d translation = transform.ExtractTranslation();
            origin = vec3(translation[0], translation[1], translation[2]);
            auto frustum = camera.GetFrustum();
            auto look_at_p = frustum.ComputeLookAtPoint();
            look_at = vec3(look_at_p[0], look_at_p[1], look_at_p[2]);
            dbg(translation, look_at_p);
        } else if (prim.IsA<pxr::UsdGeomMesh>()) {
            auto mesh = pxr::UsdGeomMesh(prim);
            pxr::VtArray<pxr::GfVec3f> points;
            pxr::VtArray<int> vertex_indices;
            pxr::VtArray<int> vertex_counts;
            mesh.GetPointsAttr().Get(&points);
            mesh.GetFaceVertexIndicesAttr().Get(&vertex_indices);
            mesh.GetFaceVertexCountsAttr().Get(&vertex_counts);

            bool valid = true;
            for (auto count : vertex_counts) {
                if (count != 3) {
                    printf("Mesh has non-triangle faces: %i\n", count);
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }

            auto transform = cache.GetLocalToWorldTransform(prim);
            pxr::GfVec3d translation = transform.ExtractTranslation();
            dbg(translation);
        }
    }

    dbg(spheres[0].center.y);
    dbg(spheres.size());

    auto scene = Scene{.spheres = std::move(spheres), .environment_map = Image("scenes/san.exr")};

    std::random_device rand_dev;
    std::default_random_engine rng(rand_dev());
    std::uniform_real_distribution<float> float_dist(0.0, 1.0);

    TevConnection connection = TevConnection();

    const std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    const uint32_t width = 960;
    const uint32_t height = 540;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    auto green = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.2, 0.9, 0.2)}}};
    auto blue = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.3, 0.3, 0.8)}}};

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

    const mat4 view = lookAt(origin, look_at, vec3(0.0, 1.0, 0.0));
    const mat4 proj = perspective(radians(fov), float(width) / float(height), 0.001f, 10000.0f);

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

                if (auto intersection_ = scene.find_intersection(ray, false)) {
                    auto [intersection, bsdf] = intersection_.value();

                    switch (bsdf.tag) {
                    case Bsdf::Tag::Diffuse: {
                        auto params = bsdf.params.diffuse;

                        auto [value, direction] =
                            scene.environment_map.sample_pdf(vec2(float_dist(rng), float_dist(rng)));

                        auto shadow_ray = Ray{intersection.position + direction * vec3(0.0001), direction, 10000.0f};

                        if (!scene.find_intersection(shadow_ray, true)) {
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

                        if (!scene.find_intersection(reflected_ray, true)) {
                            colour = scene.environment_map.sample_env_map(reflected_ray.d) * params.colour;
                        }

                        break;
                    }
                    }
                } else {
                    colour = scene.environment_map.sample_env_map(ray.d);
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
