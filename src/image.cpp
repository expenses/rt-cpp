#include "debugging.hpp"
#include "sampling.hpp"
#include "util.hpp"
#include <fstream>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <vector>

#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfInputFile.h>
#include <ImfMatrixAttribute.h>
#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>

using namespace glm;

struct Image {
    Image(const std::string &filename) {
        Imf::InputFile file(filename.data());
        const auto header = file.header();
        const auto data_window = header.dataWindow();
        const auto dimensions = data_window.max;
        width = dimensions.x + 1;
        height = dimensions.y + 1;
        data.resize(width * height * 3);
        Imf::FrameBuffer framebuffer;
        framebuffer.insert("R", Imf::Slice(Imf::FLOAT, (char *)&data[0], sizeof(float) * 3, sizeof(float) * width * 3));
        framebuffer.insert("G", Imf::Slice(Imf::FLOAT, (char *)&data[1], sizeof(float) * 3, sizeof(float) * width * 3));
        framebuffer.insert("B", Imf::Slice(Imf::FLOAT, (char *)&data[2], sizeof(float) * 3, sizeof(float) * width * 3));
        file.setFrameBuffer(framebuffer);
        file.readPixels(0, height - 1);

        std::vector<float> luminance_values(width * height);

        float theta_scale = 1.0f / float(height - 1) * M_PI;

        for (int y = 0; y < height; y++) {
            // A scaling factor, as the rows near the poles are less important
            // than
            //
            float sin_theta = sin(y * theta_scale);

            for (int x = 0; x < width; x++) {
                uint64_t offset = (y * width + x) * 3;

                auto rgb = vec3(data[offset], data[offset + 1], data[offset + 2]);

                luminance_values[y * width + x] = luminance(rgb) * sin_theta;
            }
        }

        distribution = Distribution2D(luminance_values, width, height);
    }

    const vec3 sample(vec2 uv) {
        if (uv.x < 0.0) {
            uv.x += 1.0;
        }

        auto x = uint32_t((uv.x * width));
        auto y = uint32_t((uv.y * height));
        x = min(x, width - 1);
        y = min(y, height - 1);

        auto offset = (y * width + x) * 3;

        return vec3(data[offset], data[offset + 1], data[offset + 2]);
    }

    const vec3 sample_env_map(vec3 dir) {
        return sample(vec2(atan2(-dir.z, -dir.x) / (2.0 * M_PI), acos(dir.y) / M_PI));
    }

    const std::pair<vec3, vec3> sample_pdf(vec2 rng) {
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

    uint32_t width;
    uint32_t height;
    std::vector<float> data;
    Distribution2D distribution;
};
