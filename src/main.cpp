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

vec3 pxr_to_vec3(pxr::GfVec3f vec) {
    return vec3(vec[0], vec[1], vec[2]);
}

pxr::GfVec3f vec3_to_pxr(vec3 vec) {
    return pxr::GfVec3f(vec.x, vec.y, vec.z);
}

float UniformFloat(pcg32_random_t &rng) {
    return pcg32_random_r(&rng) * 0x1p-32f;
}

int main(int argc, char *argv[]) {
    OIIO::TextureSystem &tex_sys = *OIIO::TextureSystem::create();
    // default: 1024
    tex_sys.attribute("max_memory_MB", 8096.0f);
    // default: 100
    tex_sys.attribute("max_open_files", 256);

    auto texture_path = OIIO::ustring("scenes/textures/7268504077753552595.jpg");

    RTCDevice device = rtcNewDevice(NULL);
    RTCScene rtscene = rtcNewScene(device);

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(argv[1]);
    pxr::UsdGeomXformCache cache;

    std::vector<Sphere> spheres;
    std::vector<Bsdf> sphere_bsdfs;
    std::vector<Mesh> meshes;
    auto origin = vec3(-10, 0, -10);
    mat4 view = glm::lookAt(origin, vec3(0), vec3(0, 1, 0));
    float fov = 59.0;
    std::string environment_map_path = "scenes/noon_grass_4k.exr";

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
                col = pxr_to_vec3(display_colours[0]);
            }

            float radiance;
            prim.GetAttribute(pxr::TfToken("radiance")).Get(&radiance);
            dbg(radiance);

            Bsdf bsdf;

            if (radiance > 0.0f) {
                bsdf = Bsdf{.tag = Bsdf::Tag::Emissive, .params = Bsdf::Params{.emissive = {.radiance = radiance}}};
            } else {
                bsdf = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = col}}};
            }

            spheres.push_back(Sphere{.center = pxr_to_vec3(pxr::GfVec3f(translation)), .radius = float(radius)});
            sphere_bsdfs.push_back(bsdf);
        } else if (prim.IsA<pxr::UsdGeomCamera>()) {
            pxr::GfCamera camera = pxr::UsdGeomCamera(prim).GetCamera(pxr::UsdTimeCode());
            fov = camera.GetFieldOfView(pxr::GfCamera::FOVDirection::FOVHorizontal);
            auto transform = camera.GetTransform();
            origin = pxr_to_vec3(pxr::GfVec3f(transform.ExtractTranslation())) + vec3(0, 5, 0);
            auto frustum = camera.GetFrustum();
            auto look_at = pxr_to_vec3(pxr::GfVec3f(frustum.ComputeLookAtPoint()));
            auto up = pxr_to_vec3(pxr::GfVec3f(frustum.ComputeUpVector()));
            view = glm::lookAt(origin, look_at, up);
            dbg(vec3_to_pxr(origin), vec3_to_pxr(look_at), vec3_to_pxr(up));
        } else if (prim.IsA<pxr::UsdGeomMesh>()) {
            auto context = MeshContext{.embree_device = device, .embree_scene = rtscene, .xform_cache = cache};

            if (auto mesh = load_mesh(prim, context)) {
                meshes.push_back(mesh.value());
            }
        } else if (prim.IsA<pxr::UsdLuxDomeLight>()) {
            auto dome_light = pxr::UsdLuxDomeLight(prim);
            pxr::SdfAssetPath path;
            dome_light.GetTextureFileAttr().Get(&path);
            // environment_map_path = path.GetResolvedPath();
        }
    }

    dbg(spheres.size());

    auto sphere_geometry = RTC_INVALID_GEOMETRY_ID;

    if (spheres.size() > 0) {
        RTCGeometry sphere_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

        Sphere *vb = (Sphere *)rtcSetNewGeometryBuffer(sphere_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4,
                                                       sizeof(Sphere), spheres.size());

        memcpy(vb, spheres.data(), spheres.size() * sizeof(Sphere));

        rtcCommitGeometry(sphere_geom);
        sphere_geometry = rtcAttachGeometry(rtscene, sphere_geom);
        rtcReleaseGeometry(sphere_geom);
    }
    rtcCommitScene(rtscene);

    auto scene =
        Scene{.sphere_bsdfs = std::move(sphere_bsdfs), .environment_map = Image(environment_map_path, tex_sys)};

    pcg32_random_t rng;

    TevConnection connection = TevConnection();

    const std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    const uint32_t width = 1920/4;
    const uint32_t height = 1080/4;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    const mat4 view_inverse = inverse(view);
    const mat4 proj_inverse = inverse(perspective(radians(fov), float(width) / float(height), 0.001f, 10000.0f));

    auto start = std::chrono::steady_clock::now();
    bool never_updated = true;

    OIIO::TextureOpt opt;
    opt.swrap = OIIO::TextureOpt::Wrap::WrapPeriodic;
    opt.twrap = OIIO::TextureOpt::Wrap::WrapPeriodic;

    int num_shadow_rays = 10;

    while (true) {
        tbb::parallel_for(tbb::blocked_range<int>(0, height), [&](tbb::blocked_range<int> y_range) {
            RTCIntersectContext context;
            context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
            rtcInitIntersectContext(&context);

            for (int y = y_range.begin(); y < y_range.end(); y++) {
                for (int x = 0; x < width; x++) {
                    const float u = (float(x) + UniformFloat(rng)) / width;
                    const float v = (float(y) + UniformFloat(rng)) / height;

                    auto target = proj_inverse * vec4(u * 2.0 - 1.0, v * -2.0 + 1.0, 1.0, 1.0);
                    auto local_direction = normalize(vec3(target));
                    auto direction = normalize(vec3(view_inverse * vec4(local_direction, 0.0)));

                    auto ray = Ray{origin, direction, 10000.0f};

                    auto rayhit = ray.as_embree();

                    rtcIntersect1(rtscene, &context, &rayhit);
                    vec3 colour = vec3(0.0);

                    ray = Ray(rayhit);

                    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                        vec3 normal = vec3(0, 1, 0);
                        Bsdf bsdf;
                        vec2 uv = vec2(0);
                        if (rayhit.hit.geomID == sphere_geometry) {
                            normal = normalize(vec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z));
                            bsdf = scene.sphere_bsdfs[rayhit.hit.primID];
                            uv = vec2(rayhit.hit.u, rayhit.hit.v);
                        } else {
                            auto inst_id = rayhit.hit.instID[0];

                            Mesh &mesh = meshes[inst_id];

                            vec3 colour = vec3(1.0, 1.0, 1.0);

                            rtcInterpolate0(mesh.geometry_handle, rayhit.hit.primID, rayhit.hit.u, rayhit.hit.v,
                                            RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
                            rtcInterpolate0(mesh.geometry_handle, rayhit.hit.primID, rayhit.hit.u, rayhit.hit.v,
                                            RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &uv.x, 2);

                            normal = pxr_to_vec3(mesh.normal_matrix * vec3_to_pxr(normalize(normal)));

                            if (mesh.materials.size() > 0) {
                                auto material_index = mesh.material_indices[rayhit.hit.primID];
                                auto material = mesh.materials.at(material_index);
                                if (material.path) {
                                    tex_sys.texture(material.path.value(), opt, uv.x, uv.y, 0.0, 0.0, 0.0, 0.0, 3,
                                                    (float *)&colour, nullptr, nullptr);
                                }
                            }

                            bsdf =
                                Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = colour}}};
                        }

                        switch (bsdf.tag) {
                        case Bsdf::Tag::Diffuse: {
                            auto pos = ray.o + ray.d * ray.max_t;
                            auto params = bsdf.params.diffuse;

                            for (int i = 0; i < num_shadow_rays; i++) {
                                auto cosine_dir = CosineSampleHemisphere(vec2(UniformFloat(rng), UniformFloat(rng)));
                                auto direction = rotation_matrix(normal) * cosine_dir;

                                /*auto [pdf, sample_uv, direction] =
                                    scene.environment_map.sample_pdf(vec2(UniformFloat(rng), UniformFloat(rng)));*/

                                //auto l_dot_n = std::max(dot(direction, normal), 0.0f);

                                /*if (l_dot_n == 0.0f) {
                                    continue;
                                }*/

                                
                                auto ray =
                                    Ray{pos + direction * vec3(0.0001), direction, 10000.0f}.as_embree();

                                rtcIntersect1(rtscene, &context, &ray);
                                
                                if (ray.hit.geomID == sphere_geometry) {
                                    auto bsdf = scene.sphere_bsdfs[ray.hit.primID];
                                    //if (bsdf.tag == Bsdf::Tag::Emissive) {
                                        colour += vec3(bsdf.params.emissive.radiance) * vec3(1.0 / M_PI);
                                    //}
                                }


                                /*
                                rtcOccluded1(rtscene, &context, &shadow_ray);

                                if (shadow_ray.tfar >= 0.0) {
                                    auto value = scene.environment_map.sample(sample_uv) / pdf;
                                    colour += value * params.colour * vec3(l_dot_n / M_PI);
                                }
                                */

                                /*if (std::isnan(colour.x) || std::isnan(colour.y) || std::isnan(colour.z)) {
                                    dbg("nan!", vec3_to_pxr(colour), vec3_to_pxr(normal), pdf,
                                        vec3_to_pxr(params.colour), l_dot_n);
                                }*/
                            }

                            colour /= num_shadow_rays;

                            break;
                        }
                        case Bsdf::Tag::Emissive: {
                            colour += bsdf.params.emissive.radiance;
                        }
                        }

                    } else {
                        // colour = scene.environment_map.sample_env_map(ray.d);
                    }

                    auto offset = (y * width + x) * 3;
                    accum.data[offset + 0] += colour.x;
                    accum.data[offset + 1] += colour.y;
                    accum.data[offset + 2] += colour.z;
                }
            }
        });

        accum.num_samples ++;

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
