#include "sampling.hpp"
#include "util.hpp"
#include <fstream>
#include <glm/geometric.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <vector>

using namespace glm;

struct PdfSample {
    vec3 value;
    vec3 direction;
};

struct Fexr {
    Fexr(const std::string &filename) {
        std::ifstream infile(filename);

        char header[4];

        infile.read((char *)&header, sizeof(header));

        assert(header[0] == 'F' && header[1] == 'E' && header[2] == 'X' &&
               header[3] == 'R');

        infile.read((char *)&width, sizeof(width));
        infile.read((char *)&height, sizeof(height));

        data.resize(width * height * 3);
        infile.read((char *)data.data(), width * height * 3 * 4);

        std::vector<float> luminance_values;
        luminance_values.resize(width * height);

        float theta_scale = 1.0 / (height - 1) * M_PI;

        for (int y = 0; y < height; y++) {
            // A scaling factor, as the rows near the poles are less important
            // than
            //
            float sin_theta = sin(y * theta_scale);

            for (int x = 0; x < width; x++) {
                uint64_t offset = (y * width + x) * 3;

                auto rgb =
                    vec3(data[offset], data[offset + 1], data[offset + 2]);

                luminance_values[y * width + x] = luminance(rgb) * sin_theta;
            }
        }

        distribution = Distribution2D(luminance_values, width, height);
    }

    vec3 sample(vec2 uv) {
        uint32_t x = (uv.x * (width - 1));
        uint32_t y = (uv.y * (height - 1));

        uint64_t offset = (y * width + x) * 3;

        return vec3(data[offset], data[offset + 1], data[offset + 2]);
    }

    vec3 sample_env_map(vec3 dir) {
        return sample(
            vec2(atan2(-dir.z, -dir.x) / (2.0 * M_PI), acos(dir.y) / M_PI));
    }

    PdfSample sample_pdf(vec2 rng) {
        PdfSample pdf_sample;

        auto dist_sample = distribution.sample(rng);

        auto theta = dist_sample.point.y * M_PI;

        auto phi = dist_sample.point.x * M_PI * 2.0;

        float pdf = dist_sample.pdf / (2.0 * M_PI * M_PI * sin(theta));

        pdf_sample.value = sample(dist_sample.point) / pdf;

        pdf_sample.direction =
            vec3(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
        pdf_sample.direction.z = -pdf_sample.direction.z;
        pdf_sample.direction.x = -pdf_sample.direction.x;

        return pdf_sample;
    }

    uint32_t width;
    uint32_t height;
    std::vector<float> data;
    Distribution2D distribution;
};
