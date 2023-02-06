#include <span>

inline void upload_vec3_geometry(RTCGeometry &geometry, std::span<pxr::GfVec3f> points, RTCBufferType type, uint slot = 0) {
    auto buffer = rtcSetNewGeometryBuffer(geometry, type, slot, RTC_FORMAT_FLOAT3, sizeof(pxr::GfVec3f), points.size());
    memcpy(buffer, points.data(), points.size() * sizeof(pxr::GfVec3f));
}

inline void upload_vec2_geometry(RTCGeometry &geometry, std::span<pxr::GfVec2f> points, RTCBufferType type, uint slot = 0) {
    auto buffer = rtcSetNewGeometryBuffer(geometry, type, slot, RTC_FORMAT_FLOAT2, sizeof(pxr::GfVec2f), points.size());
    memcpy(buffer, points.data(), points.size() * sizeof(pxr::GfVec2f));
}

inline void upload_index_buffer(RTCGeometry &geometry, std::span<int> indices) {
    auto buffer = rtcSetNewGeometryBuffer(geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(int) * 3,
                                          indices.size() / 3);
    memcpy(buffer, indices.data(), indices.size() * sizeof(int));
}
