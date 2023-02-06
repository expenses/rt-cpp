#pragma once
#include "bsdf.hpp"
#include "geo.hpp"
#include "image.hpp"

struct SphereEmitter {
    Sphere sphere;
    size_t index;
};

struct Emitters {
    std::optional<Image> environment_map;
    std::vector<SphereEmitter> sphere_emitters;
};

struct Intersection {
    vec3 normal;
    Bsdf bsdf;
};

struct Scene {
    std::vector<Bsdf> sphere_bsdfs;
    std::vector<Mesh> meshes;
    Emitters emitters;
    uint sphere_geometry_ref;
    RTCScene embree_scene;

    vec3 sample_emitter(RayHit rayhit) {
        vec3 sample = vec3(0.0);

        if (sphere_bsdfs.size() > 0 && rayhit.hit.geomID == sphere_geometry_ref) {
            auto bsdf = sphere_bsdfs[rayhit.hit.primID];
            if (bsdf.tag == Bsdf::Tag::Emissive) {
                sample = bsdf.params.emissive.radiance;
            }
        } else if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            if (emitters.environment_map) {
                sample = emitters.environment_map.value().sample_env_map(rayhit.ray.direction);
            }
        }

        return sample;
    }

    std::tuple<vec3, vec3> sample_emitter_pdf(vec3 pos, vec2 rng) {
        if (emitters.environment_map) {
            Image& environment_map = emitters.environment_map.value();

            auto [pdf, sample_point, direction] = environment_map.sample_pdf(rng);

            auto test_ray = Ray {
                .origin = pos + direction * vec3(0.0001f),
                .direction = direction,
                .t_far = 10000.0f
            };

            RTCIntersectContext context;
            rtcInitIntersectContext(&context);
            rtcOccluded1(embree_scene, &context, reinterpret_cast<RTCRay*>(&test_ray));

            if (std::isinf(test_ray.t_far)) {
                return std::make_tuple(vec3(0), direction);
            } else {
                return std::make_tuple(environment_map.sample(sample_point) / pdf, direction);
            }
        } else {
            return std::make_tuple(vec3(0), vec3(0,1,0));
        }
    }

    Intersection get_intersection_for_rayhit(RayHit rayhit, OIIO::TextureSystem &tex_sys) {
        OIIO::TextureOpt opt;
        opt.swrap = OIIO::TextureOpt::Wrap::WrapPeriodic;
        opt.twrap = OIIO::TextureOpt::Wrap::WrapPeriodic;

        if (rayhit.hit.geomID == sphere_geometry_ref) {
            auto normal = normalize(rayhit.hit.normal);
            auto bsdf = sphere_bsdfs[rayhit.hit.primID];
            return Intersection {
                .normal = normal, .bsdf = bsdf
            };
        } else {
            Mesh &mesh = meshes.at(rayhit.hit.instID[0]);

            vec3 colour = vec3(1.0, 1.0, 1.0);
            vec2 uv = vec2(0);
            vec3 normal = vec3(0);

            rtcInterpolate0(mesh.geometry_handle, rayhit.hit.primID, rayhit.hit.uv.x, rayhit.hit.uv.y,
                            RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
            rtcInterpolate0(mesh.geometry_handle, rayhit.hit.primID, rayhit.hit.uv.x, rayhit.hit.uv.y,
                            RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &uv.x, 2);

            normal = mesh.normal_matrix * normalize(normal);

            if (mesh.materials.size() > 0) {
                auto material_index = mesh.material_indices.at(rayhit.hit.primID);
                auto material = mesh.materials.at(size_t(material_index));
                if (material.path) {
                    tex_sys.texture(material.path.value(), opt, uv.x, uv.y, 0.0, 0.0, 0.0, 0.0, 3,
                                    reinterpret_cast<float *>(&colour), nullptr, nullptr);
                }
            }

            Bsdf bsdf =
            Bsdf{.tag = Bsdf::Tag::Diffuse, .params = Bsdf::Params{.diffuse = {.colour = colour}}};

            return Intersection {
                .normal = normal, .bsdf = bsdf
            };
        }
    }
};
