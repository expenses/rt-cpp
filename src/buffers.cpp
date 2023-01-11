#include <cstdint>
#include <vector>

struct OutputBuffer {
    std::vector<float> data;

    OutputBuffer(uint32_t width, uint32_t height) {
        data = std::vector<float>(width * height * 3, 0.0f);
    }
};

struct AccumulationBuffer {
    std::vector<float> data;
    uint64_t num_samples;

    AccumulationBuffer(uint32_t width, uint32_t height) {
        data = std::vector<float>(width * height * 3, 0.0f);
        num_samples = 0;
    }

    void update_output(OutputBuffer &output) {
        for (int i = 0; i < data.size(); i++) {
            output.data[i] = data[i] / num_samples;
        }
    }
};
