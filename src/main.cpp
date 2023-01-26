#include <cmath>
#include <iostream>

#include "external.hpp"

#include "buffers.hpp"
#include "debugging.hpp"
#include "geo.hpp"
#include "image.hpp"
#include "pcg.h"
#include "sampling.hpp"
#include "scene.hpp"
#include "tev.hpp"
#include "util.hpp"
#include <pxr/pxr.h>

using namespace glm;

pxr::GfVec3f pxr_vec3(vec3 vec) {
    return pxr::GfVec3f(vec.x, vec.y, vec.z);
}

float UniformFloat(pcg32_random_t &rng) {
    return pcg32_random_r(&rng) * 0x1p-32f;
}

void upload_vec3_geometry(RTCGeometry &geometry, std::span<pxr::GfVec3f> points, RTCBufferType type, int slot = 0) {
    auto buffer = rtcSetNewGeometryBuffer(geometry, type, slot, RTC_FORMAT_FLOAT3, sizeof(pxr::GfVec3f), points.size());
    memcpy(buffer, points.data(), points.size() * sizeof(pxr::GfVec3f));
}

void upload_vec2_geometry(RTCGeometry &geometry, std::span<pxr::GfVec2f> points, RTCBufferType type, int slot = 0) {
    auto buffer = rtcSetNewGeometryBuffer(geometry, type, slot, RTC_FORMAT_FLOAT2, sizeof(pxr::GfVec2f), points.size());
    memcpy(buffer, points.data(), points.size() * sizeof(pxr::GfVec2f));
}

std::string resolve_diffuse_material_texture_path(pxr::UsdShadeMaterial& material) {
    auto source = material.ComputeSurfaceSource();
    auto diffuse = source.GetInput(pxr::TfToken("diffuseColor"));
    source = diffuse.GetConnectedSources()[0].source;
    auto texture = source.GetInput(pxr::TfToken("file"));
    pxr::SdfAssetPath path;
    texture.Get(&path);
    return path.GetResolvedPath();
}

struct Mesh {
    pxr::GfMatrix3f normal_matrix;
    RTCScene scene;
    std::vector<int> material_indices;
    std::vector<Material> materials;
};

int main(int argc, char *argv[]) {
    auto tex_sys = OIIO::TextureSystem::create();

    auto texture_path = OIIO::ustring("scenes/textures/7268504077753552595.jpg");

    RTCDevice device = rtcNewDevice(NULL);
    RTCScene rtscene = rtcNewScene(device);

    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(argv[1]);
    pxr::UsdGeomXformCache cache;

    std::vector<Sphere> spheres;
    std::vector<Bsdf> sphere_bsdfs;
    std::vector<Mesh> meshes;
    vec3 origin = vec3(-1.0, 0.0, 0.0);
    vec3 look_at = vec3(0.0, 0.0, 0.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
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
                col = vec3(display_colours[0][0], display_colours[0][1], display_colours[0][2]);
            }

            auto diffuse = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = col}}};

            spheres.push_back(
                Sphere{.center = vec3(translation[0], translation[1], translation[2]), .radius = float(radius)});
            sphere_bsdfs.push_back(diffuse);
        } else if (prim.IsA<pxr::UsdGeomCamera>()) {
            pxr::GfCamera camera = pxr::UsdGeomCamera(prim).GetCamera(pxr::UsdTimeCode());
            fov = camera.GetFieldOfView(pxr::GfCamera::FOVDirection::FOVHorizontal);
            auto transform = camera.GetTransform();
            pxr::GfVec3d translation = transform.ExtractTranslation();
            origin = vec3(translation[0], translation[1], translation[2]); // + vec3(0,5,0);
            auto frustum = camera.GetFrustum();
            auto look_at_p = frustum.ComputeLookAtPoint();
            auto up_p = frustum.ComputeUpVector();
            look_at = vec3(look_at_p[0], look_at_p[1], look_at_p[2]);
            up = vec3(up_p[0], up_p[1], up_p[2]);
            dbg(up_p);
            dbg(translation, look_at_p);
        } else if (prim.IsA<pxr::UsdGeomMesh>()) {
            auto mesh = pxr::UsdGeomMesh(prim);
            auto primvars = pxr::UsdGeomPrimvarsAPI(prim);
            pxr::VtArray<pxr::GfVec3f> points;
            pxr::VtArray<int> vertex_indices;
            pxr::VtArray<int> vertex_counts;
            pxr::VtArray<pxr::GfVec3f> normals;
            pxr::VtArray<pxr::GfVec2f> uvs;
            mesh.GetPointsAttr().Get(&points);
            mesh.GetFaceVertexIndicesAttr().Get(&vertex_indices);
            mesh.GetFaceVertexCountsAttr().Get(&vertex_counts);
            mesh.GetNormalsAttr().Get(&normals);
            primvars.GetPrimvar(pxr::TfToken("primvars:UVMap")).Get(&uvs);

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

            // Blender exported normals don't seem to be indexed??
            std::vector<pxr::GfVec3f> indexed_normals;
            indexed_normals.resize(points.size());
            for (int i = 0; i < normals.size(); i++) {
                indexed_normals.at(vertex_indices[i]) = normals[i];
            }

            std::vector<pxr::GfVec2f> indexed_uvs;
            indexed_uvs.resize(points.size());
            for (int i = 0; i < uvs.size(); i++) {
                indexed_uvs.at(vertex_indices[i]) = uvs[i];
            }

            dbg(uvs.size(), normals.size(), indexed_normals.size(), indexed_uvs.size());

            RTCGeometry triangle_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

            rtcSetGeometryVertexAttributeCount(triangle_geom, 2);

            upload_vec3_geometry(triangle_geom, points, RTC_BUFFER_TYPE_VERTEX);
            upload_vec3_geometry(triangle_geom, indexed_normals, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0);
            upload_vec2_geometry(triangle_geom, indexed_uvs, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1);

            auto buf2 = rtcSetNewGeometryBuffer(triangle_geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
                                                sizeof(int) * 3, vertex_indices.size() / 3);
            memcpy(buf2, vertex_indices.data(), vertex_indices.size() * sizeof(int));

            std::vector<int> material_indices(vertex_indices.size(), 0);

            auto subsets = pxr::UsdGeomSubset::GetAllGeomSubsets(mesh);

            std::vector<Material> materials;

            for (int subset_id = 0; subset_id < subsets.size(); subset_id++) {
                auto material =
                    pxr::UsdShadeMaterialBindingAPI(subsets[subset_id]).ComputeBoundMaterial(pxr::UsdShadeTokens->full);
                pxr::VtArray<int> subset_indices;
                subsets[subset_id].GetIndicesAttr().Get(&subset_indices);

                for (auto i : subset_indices) {
                    material_indices[i] = subset_id;
                }

                
                auto path = resolve_diffuse_material_texture_path(material);

                materials.push_back(Material {
                    .tag = Bsdf::Tag::Diffuse,
                    .path = OIIO::ustring(path)
                });
            }

            if (subsets.size() == 0) {
                auto material =
                    pxr::UsdShadeMaterialBindingAPI(mesh).ComputeBoundMaterial(pxr::UsdShadeTokens->full);
                auto path = resolve_diffuse_material_texture_path(material);
                materials.push_back(Material {
                    .tag = Bsdf::Tag::Diffuse,
                    .path = OIIO::ustring(path)
                });
            }

            RTCScene mesh_scene = rtcNewScene(device);
            rtcCommitGeometry(triangle_geom);
            rtcAttachGeometry(mesh_scene, triangle_geom);
            rtcReleaseGeometry(triangle_geom);
            rtcCommitScene(mesh_scene);

            auto transform = pxr::GfMatrix4f(cache.GetLocalToWorldTransform(prim));

            RTCGeometry instance_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_INSTANCE);
            rtcSetGeometryInstancedScene(instance_geom, mesh_scene);
            rtcSetGeometryTransform(instance_geom, 0, RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR, transform.data());
            rtcSetGeometryTimeStepCount(instance_geom, 1);
            rtcCommitGeometry(instance_geom);
            rtcAttachGeometry(rtscene, instance_geom);
            rtcReleaseGeometry(instance_geom);

            meshes.push_back(
                Mesh{.normal_matrix = transform.GetOrthonormalized().GetTranspose().ExtractRotationMatrix(),
                     .scene = mesh_scene,
                     .material_indices = material_indices,
                     .materials = materials});
        } else if (prim.IsA<pxr::UsdLuxDomeLight>()) {
            auto dome_light = pxr::UsdLuxDomeLight(prim);
            pxr::SdfAssetPath path;
            dome_light.GetTextureFileAttr().Get(&path);
            //environment_map_path = path.GetResolvedPath();
        }
    }

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

    auto scene = Scene{.sphere_bsdfs = std::move(sphere_bsdfs),
                       .environment_map = Image(environment_map_path)};

    pcg32_random_t rng;

    TevConnection connection = TevConnection();

    const std::string image_name = "output";
    std::vector<std::string> channel_names = {"R", "G", "B"};
    std::vector<uint64_t> channel_offsets = {0, 1, 2};
    std::vector<uint64_t> channel_strides = {3, 3, 3};

    const uint32_t width = 256;
    const uint32_t height = 256;

    auto accum = AccumulationBuffer(width, height);
    auto output = OutputBuffer(width, height);

    connection.send_create(image_name, width, height, channel_names);

    auto grey = Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = vec3(0.5, 0.5, 0.5)}}};
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

    const mat4 view = lookAt(origin, look_at, up);
    const mat4 proj = perspective(radians(fov/8), float(width) / float(height), 0.001f, 10000.0f);

    const mat4 view_inverse = inverse(view);
    const mat4 proj_inverse = inverse(proj);

    auto start = std::chrono::steady_clock::now();
    bool never_updated = true;

    std::vector<RTCRayHit> rays;

    RTCIntersectContext context;
    context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
    rtcInitIntersectContext(&context);

    OIIO::TextureOpt opt;
    opt.swrap = OIIO::TextureOpt::Wrap::WrapPeriodic;
    opt.twrap = OIIO::TextureOpt::Wrap::WrapPeriodic;
    

    while (true) {
        rays.clear();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                const float u = (float(x) + UniformFloat(rng)) / width;
                const float v = (float(y) + UniformFloat(rng)) / height;

                auto target = proj_inverse * vec4(u * 2.0 - 1.0, v * -2.0 + 1.0, 1.0, 1.0);
                auto local_direction = normalize(vec3(target));
                auto direction = normalize(vec3(view_inverse * vec4(local_direction, 0.0)));

                auto ray = Ray{origin, direction, 10000.0f};

                rays.push_back(ray.as_embree());
            }
        }

        rtcIntersect1M(rtscene, &context, rays.data(), rays.size(), sizeof(RTCRayHit));

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                vec3 colour = vec3(0.0);

                auto rayhit = rays[(y * width + x)];

                auto ray = Ray(rayhit);

                if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                    vec3 normal;
                    Bsdf bsdf;
                    vec2 uv = vec2(0);
                    if (rayhit.hit.geomID == sphere_geometry) {
                        normal = normalize(vec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z));
                        bsdf = scene.sphere_bsdfs[rayhit.hit.primID];
                        uv = vec2(rayhit.hit.u, rayhit.hit.v);
                    } else {
                        auto inst_id = rayhit.hit.instID[0];

                        Mesh& mesh = meshes[inst_id];

                        vec3 colour = vec3(1.0, 0.0, 1.0);

                        rtcInterpolate0(rtcGetGeometry(mesh.scene, 0), rayhit.hit.primID, rayhit.hit.u,
                                        rayhit.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
                        rtcInterpolate0(rtcGetGeometry(mesh.scene, 0), rayhit.hit.primID, rayhit.hit.u,
                                        rayhit.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &uv.x, 2);

                        auto n = mesh.normal_matrix * pxr::GfVec3f(normal.x, normal.y, normal.z);
                        normal = vec3(n[0], n[1], n[2]);
                        
                        if (mesh.materials.size() > 0) {
                            auto material_index = mesh.material_indices[rayhit.hit.primID];
                            auto material = mesh.materials.at(material_index);
                            tex_sys->texture(material.path, opt, uv.x, uv.y, 0.0, 0.0, 0.0, 0.0, 3, (float*)&colour, nullptr, nullptr);
                        }

                        bsdf = Bsdf{.tag = Bsdf::Tag::Diffuse,
                                    .params = Bsdf::Params{.diffuse = {.colour = colour}}};
                    }

                    auto pos = ray.o + ray.d * ray.max_t;

                    auto [value, direction] =
                        scene.environment_map.sample_pdf(vec2(UniformFloat(rng), UniformFloat(rng)));

                    auto params = bsdf.params.diffuse;

                    auto shadow_ray = Ray{pos + direction * vec3(0.0001), direction, 10000.0f}.as_embree().ray;

                    rtcOccluded1(rtscene, &context, &shadow_ray);

                    auto l_dot_n = std::max(dot(direction, normal), 0.0f);

                    if (shadow_ray.tfar >= 0.0) {
                        colour = value * params.colour * vec3(l_dot_n / M_PI);
                    }

                    if (std::isnan(colour.x) || std::isnan(colour.y) || std::isnan(colour.z)) {
                        dbg("nan!", pxr_vec3(colour), pxr_vec3(value), pxr_vec3(params.colour), l_dot_n);
                    }

                    //colour = vec3(u_deriv.x, u_deriv.y, 0.0);
                    //colour = vec3(uv.x, uv.y, 0.0);
                } else {
                    colour = scene.environment_map.sample_env_map(ray.d);
                }

                /*if (auto intersection_ = scene.find_intersection(ray, false)) {
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
                }*/

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
