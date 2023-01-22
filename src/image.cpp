#include "image.hpp"
#include "debugging.hpp"
#include "external.hpp"
#include "util.hpp"
#include <vector>

using namespace glm;

Image::Image(const std::string &filename) {
    Imf::InputFile file(filename.data());
    const auto header = file.header();
    const auto data_window = header.dataWindow();
    const auto dimensions = data_window.max;
    width = dimensions.x + 1;
    height = dimensions.y + 1;
    rgb.resize(width * height);
    Imf::FrameBuffer framebuffer;
    const auto row_stride = sizeof(float) * 3 * width;
    const auto interlacing = sizeof(float) * 3;
    framebuffer.insert("R", Imf::Slice(Imf::FLOAT, (char *)rgb.data(), interlacing, row_stride));
    framebuffer.insert("G", Imf::Slice(Imf::FLOAT, (char *)rgb.data() + sizeof(float), interlacing, row_stride));
    framebuffer.insert("B", Imf::Slice(Imf::FLOAT, (char *)rgb.data() + sizeof(float) * 2, interlacing, row_stride));
    file.setFrameBuffer(framebuffer);
    file.readPixels(0, height - 1);

    std::vector<float> luminance_values(width * height);

    float theta_scale = 1.0 / float(height - 1) * M_PI;

    for (int y = 0; y < height; y++) {
        float sin_theta = sin(y * theta_scale);

        for (int x = 0; x < width; x++) {
            uint64_t offset = y * width + x;
            luminance_values[offset] = luminance(rgb[offset]) * sin_theta;
        }
    }

    distribution = Distribution2D(luminance_values, width, height);
}

const vec3 Image::sample(vec2 uv) {
    if (uv.x < 0.0) {
        uv.x += 1.0;
    }

    auto x = uint32_t((uv.x * width));
    auto y = uint32_t((uv.y * height));
    x = min(x, width - 1);
    y = min(y, height - 1);

    return rgb[y * width + x];
}

const vec3 Image::sample_env_map(vec3 dir) {
    return sample(vec2(atan2(-dir.z, -dir.x) / (2.0 * M_PI), acos(dir.y) / M_PI));
}

const std::pair<vec3, vec3> Image::sample_pdf(vec2 rng) {
    auto dist_sample = distribution.sample(rng);

    auto theta = dist_sample.point.y * M_PI;

    auto phi = dist_sample.point.x * M_PI * 2.0;

    float pdf = dist_sample.pdf / (2.0 * M_PI * M_PI * sin(theta));

    auto value = sample(dist_sample.point) / pdf;

    auto direction = vec3(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
    direction.z = -direction.z;
    direction.x = -direction.x;

    return std::make_pair(value, direction);
}
