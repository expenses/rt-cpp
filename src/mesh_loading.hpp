#include "bsdf.hpp"

struct MeshContext {
    RTCDevice &embree_device;
    RTCScene &embree_scene;
    pxr::UsdGeomXformCache &xform_cache;
};

struct Mesh {
    pxr::GfMatrix3f normal_matrix;
    RTCScene scene;
    std::vector<int> material_indices;
    std::vector<Material> materials;
    RTCGeometry geometry_handle;
};

std::optional<std::string> resolve_diffuse_material_texture_path(pxr::UsdShadeMaterial &material) {
    auto source = material.ComputeSurfaceSource();
    if (!source) {
        return std::nullopt;
    }
    auto diffuse = source.GetInput(pxr::TfToken("diffuseColor"));
    auto sources = diffuse.GetConnectedSources();

    if (sources.size() == 0) {
        return std::nullopt;
    }

    source = sources[0].source;
    auto texture = source.GetInput(pxr::TfToken("file"));
    pxr::SdfAssetPath path;
    texture.Get(&path);
    return path.GetResolvedPath();
}

std::optional<Mesh> load_mesh(pxr::UsdPrim &prim, MeshContext &context) {
    auto mesh = pxr::UsdGeomMesh(prim);

    auto primvars = pxr::UsdGeomPrimvarsAPI(mesh);
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
    auto vertex_offset = 0;
    std::vector<int> trianglulated_indices;
    for (auto count : vertex_counts) {
        if (count > 4) {
            printf("Mesh has non-triangle faces: %i\n", count);
            return std::nullopt;
        } else if (count == 3) {
            trianglulated_indices.push_back(vertex_indices[vertex_offset]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 1]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 2]);
        } else if (count == 4) {
            trianglulated_indices.push_back(vertex_indices[vertex_offset]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 1]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 2]);

            trianglulated_indices.push_back(vertex_indices[vertex_offset]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 2]);
            trianglulated_indices.push_back(vertex_indices[vertex_offset + 3]);
        }

        vertex_offset += count;
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

    /*dbg(points.size(), uvs.size(), normals.size(), indexed_normals.size(), indexed_uvs.size(),
        trianglulated_indices.size(), vertex_indices.size());*/

    RTCGeometry triangle_geom = rtcNewGeometry(context.embree_device, RTC_GEOMETRY_TYPE_TRIANGLE);

    rtcSetGeometryVertexAttributeCount(triangle_geom, 2);

    upload_vec3_geometry(triangle_geom, points, RTC_BUFFER_TYPE_VERTEX);
    upload_vec3_geometry(triangle_geom, indexed_normals, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0);
    upload_vec2_geometry(triangle_geom, indexed_uvs, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1);
    upload_index_buffer(triangle_geom, trianglulated_indices);

    std::vector<int> material_indices(trianglulated_indices.size(), 0);

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

        std::optional<OIIO::ustring> oiio_path = std::nullopt;
        if (auto path = resolve_diffuse_material_texture_path(material)) {
            oiio_path = OIIO::ustring(path.value());
        }

        materials.push_back(Material{.tag = Bsdf::Tag::Diffuse, .path = oiio_path});
    }

    if (subsets.size() == 0) {
        auto material = pxr::UsdShadeMaterialBindingAPI(mesh).ComputeBoundMaterial(pxr::UsdShadeTokens->full);

        std::optional<OIIO::ustring> oiio_path = std::nullopt;
        if (auto path = resolve_diffuse_material_texture_path(material)) {
            oiio_path = OIIO::ustring(path.value());
        }
        materials.push_back(Material{.tag = Bsdf::Tag::Diffuse, .path = oiio_path});
    }

    RTCScene mesh_scene = rtcNewScene(context.embree_device);
    rtcCommitGeometry(triangle_geom);
    rtcAttachGeometry(mesh_scene, triangle_geom);
    rtcReleaseGeometry(triangle_geom);
    rtcCommitScene(mesh_scene);

    auto transform = pxr::GfMatrix4f(context.xform_cache.GetLocalToWorldTransform(prim));

    RTCGeometry instance_geom = rtcNewGeometry(context.embree_device, RTC_GEOMETRY_TYPE_INSTANCE);
    rtcSetGeometryInstancedScene(instance_geom, mesh_scene);
    rtcSetGeometryTransform(instance_geom, 0, RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR, transform.data());
    rtcSetGeometryTimeStepCount(instance_geom, 1);
    rtcCommitGeometry(instance_geom);
    rtcAttachGeometry(context.embree_scene, instance_geom);
    rtcReleaseGeometry(instance_geom);

    auto handle = rtcGetGeometry(mesh_scene, 0);

    return Mesh{.normal_matrix = transform.GetOrthonormalized().GetTranspose().ExtractRotationMatrix(),
                .scene = mesh_scene,
                .material_indices = material_indices,
                .materials = materials,
                .geometry_handle = handle};
}
