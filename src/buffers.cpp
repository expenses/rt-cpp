#include "buffers.hpp"

OutputBuffer::OutputBuffer(uint32_t width, uint32_t height)
    : data(std::vector<float>(uint64_t(width * height) * 3, 0.0f)) {
}

AccumulationBuffer::AccumulationBuffer(uint32_t width, uint32_t height)
    : data(std::vector<float>(uint64_t(width * height) * 3, 0.0f)) {
}

void AccumulationBuffer::update_output(OutputBuffer &output) {
    for (int i = 0; i < data.size(); i++) {
        output.data[i] = data[i] / float(num_samples);
    }
}
