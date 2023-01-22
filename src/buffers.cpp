#include "buffers.h"

OutputBuffer::OutputBuffer(uint32_t width, uint32_t height) {
    data = std::vector<float>(width * height * 3, 0.0f);
}

AccumulationBuffer::AccumulationBuffer(uint32_t width, uint32_t height) {
    data = std::vector<float>(width * height * 3, 0.0f);
    num_samples = 0;
}

void AccumulationBuffer::update_output(OutputBuffer &output) {
    for (int i = 0; i < data.size(); i++) {
        output.data[i] = data[i] / num_samples;
    }
}
