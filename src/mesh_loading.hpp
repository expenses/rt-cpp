#include "bsdf.hpp"

struct MeshContext {
    RTCDevice &embree_device;
    RTCScene &embree_scene;
    pxr::UsdGeomXformCache &xform_cache;
};

struct Mesh {
    mat3 normal_matrix;
    RTCScene scene;
    std::vector<int> material_indices;
    std::vector<Material> materials;
    RTCGeometry geometry_handle;
};

std::optional<Mesh> load_mesh(pxr::UsdPrim &prim, MeshContext &context);
