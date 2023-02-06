#include <cmath>
#include <iostream>

#include "external.hpp"

#include "buffers.hpp"
#include "debugging.hpp"
#include "embree_helpers.hpp"
#include "geo.hpp"
#include "image.hpp"
#include "mesh_loading.hpp"
#include "pcg.h"
#include "sampling.hpp"
#include "scene.hpp"
#include "tev.hpp"
#include "util.hpp"

using namespace glm;

struct CustomAttributes {
    const pxr::TfToken RADIANCE = pxr::TfToken("custom_radiance");
    const pxr::TfToken MIRROR_COLOUR = pxr::TfToken("custom_mirror_colour");

    vec3 radiance = vec3(0);
    vec3 mirror_colour = vec3(0);

    CustomAttributes(pxr::UsdPrim &prim) {
        pxr::GfVec3f radiance_pxr = pxr::GfVec3f(0);
        prim.GetAttribute(RADIANCE).Get(&radiance_pxr);
        radiance = pxr_to_vec3(radiance_pxr);
        
        pxr::GfVec3f mirror_colour_pxr = pxr::GfVec3f(0);
        prim.GetAttribute(MIRROR_COLOUR).Get(&mirror_colour_pxr);
        dbg(mirror_colour_pxr);
        mirror_colour = pxr_to_vec3(mirror_colour_pxr);
        
    }
};

static float UniformFloat(pcg32_random_t &rng) {
    return float(pcg32_random_r(&rng)) * 0x1p-32f;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Expected 1 param, got %i\n", argc - 1);
        return -1;
    }

    OIIO::TextureSystem &tex_sys = *OIIO::TextureSystem::create();
    // default: 1024
    tex_sys.attribute("max_memory_MB", 8096.0f);
    // default: 100
    tex_sys.attribute("max_open_files", 256);

    auto texture_path = OIIO::ustring("scenes/textures/7268504077753552595.jpg");

    RTCDevice device = rtcNewDevice(nullptr);

    Scene scene;
    scene.embree_scene = rtcNewScene(device);

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(argv[1]);
    pxr::UsdGeomXformCache cache;

    std::vector<Sphere> spheres;
    auto origin = vec3(-10, 0, -10);
    mat4 view = glm::lookAt(origin, vec3(0), vec3(0, 1, 0));
    float fov = 59.0;

    for (pxr::UsdPrim prim : pxr::UsdPrimRange::Stage(stage)) {
        if (prim.IsA<pxr::UsdGeomSphere>()) {
            auto usd_sphere = pxr::UsdGeomSphere(prim);

            double radius;
            usd_sphere.GetRadiusAttr().Get(&radius);

            auto transform = cache.GetLocalToWorldTransform(prim);
            pxr::GfVec3d translation = transform.ExtractTranslation();

            pxr::VtArray<pxr::GfVec3f> display_colours;
            prim.GetAttribute(pxr::UsdGeomTokensType().primvarsDisplayColor).Get(&display_colours);
            vec3 col = vec3(0.5, 0.5, 0.5);
            if (display_colours.size() > 0) {
                col = pxr_to_vec3(display_colours[0]);
            }

            CustomAttributes custom_attributes = CustomAttributes(prim);

            Bsdf bsdf;
            Sphere sphere = Sphere{.center = pxr_to_vec3(pxr::GfVec3f(translation)), .radius = float(radius)};
            spheres.push_back(sphere);

            if (custom_attributes.radiance != vec3(0.0)) {
                bsdf = Bsdf{.tag = Bsdf::Tag::Emissive,
                            .params = Bsdf::Params{.emissive = {.radiance = custom_attributes.radiance}}};
            } else if (custom_attributes.mirror_colour != vec3(0.0)) {
                bsdf = Bsdf{.tag = Bsdf::Tag::Mirror, .params = Bsdf::Params{.mirror = {.colour = custom_attributes.mirror_colour}}};
            } else {
                bsdf = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = col}}};
            }

            scene.sphere_bsdfs.push_back(bsdf);

            if (bsdf.tag == Bsdf::Tag::Emissive) {
                scene.emitters.sphere_emitters.push_back(SphereEmitter {
                    .sphere = sphere,
                    .index = scene.sphere_bsdfs.size() - 1
                });
            }
        } else if (prim.IsA<pxr::UsdGeomCamera>()) {
            pxr::GfCamera camera = pxr::UsdGeomCamera(prim).GetCamera(pxr::UsdTimeCode());
            fov = camera.GetFieldOfView(pxr::GfCamera::FOVDirection::FOVHorizontal) / 3.0;
            auto transform = camera.GetTransform();
            origin = pxr_to_vec3(pxr::GfVec3f(transform.ExtractTranslation()));
            auto frustum = camera.GetFrustum();
            auto look_at = pxr_to_vec3(pxr::GfVec3f(frustum.ComputeLookAtPoint())) + vec3(0,-3.0,0);
            auto up = pxr_to_vec3(pxr::GfVec3f(frustum.ComputeUpVector()));
            view = glm::lookAt(origin, look_at, up);
        } else if (prim.IsA<pxr::UsdGeomMesh>()) {
            auto context = MeshContext{.embree_device = device, .embree_scene = scene.embree_scene, .xform_cache = cache};

            if (auto mesh = load_mesh(prim, context)) {
                scene.meshes.push_back(mesh.value());
            }
        } else if (prim.IsA<pxr::UsdLuxDomeLight>()) {
            auto dome_light = pxr::UsdLuxDomeLight(prim);
            pxr::SdfAssetPath path;
            dome_light.GetTextureFileAttr().Get(&path);
            scene.emitters.environment_map.emplace(Image(path.GetResolvedPath(), tex_sys));
        } else {
            auto string = prim.GetTypeName();
            if (string == "Xform" || string == "Shader" || string == "Material" || string == "Scope" ||
                string == "GeomSubset" || string == "") {
                // Ignored.
            } else {
                dbg(prim.GetName(), prim.GetTypeName());
            }
        }
    }

    dbg(spheres.size());

    scene.sphere_geometry_ref = RTC_INVALID_GEOMETRY_ID;

    if (spheres.size() > 0) {
        RTCGeometry sphere_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

        Sphere *vb = reinterpret_cast<Sphere *>(rtcSetNewGeometryBuffer(
            sphere_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, sizeof(Sphere), spheres.size()));

        memcpy(vb, spheres.data(), spheres.size() * sizeof(Sphere));

        rtcCommitGeometry(sphere_geom);
        scene.sphere_geometry_ref = rtcAttachGeometry(scene.embree_scene, sphere_geom);
        rtcReleaseGeometry(sphere_geom);
    }
    rtcCommitScene(scene.embree_scene);

    pcg32_random_t rng;

    TevConnection connection = TevConnection();

    const std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    const uint32_t width = 1920 / 4;
    const uint32_t height = 1080 / 4;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    const mat4 view_inverse = inverse(view);
    const mat4 proj_inverse = inverse(perspective(radians(fov), float(width) / float(height), 0.001f, 10000.0f));

    auto start = std::chrono::steady_clock::now();
    bool never_updated = true;

    int num_shadow_rays = 10;

    while (true) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size_t(height)), [&](tbb::blocked_range<size_t> y_range) {
            RTCIntersectContext context;
            context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
            rtcInitIntersectContext(&context);

            //for (size_t y = 0; y < height; y++) {
            for (size_t y = y_range.begin(); y < y_range.end(); y++) {
                for (size_t x = 0; x < width; x++) {
                    const float u = (float(x) + UniformFloat(rng)) / width;
                    const float v = (float(y) + UniformFloat(rng)) / height;

                    auto target = proj_inverse * vec4(u * 2.0f - 1.0f, v * -2.0f + 1.0f, 1.0f, 1.0f);
                    auto local_direction = normalize(vec3(target));
                    auto direction = normalize(vec3(view_inverse * vec4(local_direction, 0.0f)));

                    auto rayhit = RayHit {
                        .ray = Ray {
                            .origin = origin,
                            .direction = direction,
                            .t_far = 10000.0f
                        },
                        .hit = Hit {
                            .geomID = RTC_INVALID_GEOMETRY_ID
                        }
                    };
                    
                    vec3 colour = vec3(0.0);

                    vec3 absortion = vec3(1.0);

                    bool done_looping = false;

                    for (int path_i = 0; path_i < 6; path_i++) {
                        if (done_looping) {
                            break;
                        }

                        rtcIntersect1(scene.embree_scene, &context, reinterpret_cast<RTCRayHit*>(&rayhit));

                        auto ray = rayhit.ray;
                        auto pos = ray.at(ray.t_far);
                    
                        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                            auto [normal, bsdf] = scene.get_intersection_for_rayhit(rayhit, tex_sys);

                            switch (bsdf.tag) {
                            case Bsdf::Tag::Diffuse: {
                                for (int i = 0; i < num_shadow_rays; i++) {
                                    auto [sample, direction] = scene.sample_emitter_pdf(pos, vec2(UniformFloat(rng), UniformFloat(rng)));
                                    colour += sample * eval_diffuse(bsdf.params.diffuse, direction, normal);
                                }

                                colour /= num_shadow_rays;
                                colour *= absortion;

                                done_looping = true;

                                break;
                            }
                            case Bsdf::Tag::Mirror: {
                                auto reflected = reflect(ray.direction, normal);

                                rayhit = RayHit {
                                    .ray = Ray {
                                        .origin = ray.at(ray.t_far) + reflected * 0.0001f,
                                        .direction = reflected,
                                        .t_far = 10000.0f
                                    },
                                    .hit = Hit {
                                        .geomID = RTC_INVALID_GEOMETRY_ID
                                    }
                                };

                                absortion *= bsdf.params.mirror.colour;

                                break;
                            }
                            case Bsdf::Tag::Emissive: {
                                colour += bsdf.params.emissive.radiance * absortion;

                                done_looping = true;

                                break;
                            }
                            }

                        } else {
                            if (scene.emitters.environment_map) {
                                colour = scene.emitters.environment_map.value().sample_env_map(rayhit.ray.direction) * absortion;
                            }
                        }
                    }

                    size_t offset = (y * width + x) * 3;
                    accum.data[offset + 0] += colour.x;
                    accum.data[offset + 1] += colour.y;
                    accum.data[offset + 2] += colour.z;
                }
            }
        });

        accum.num_samples++;

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
