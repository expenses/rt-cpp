#pragma once
#include <cstdint>
#include <vector>

struct OutputBuffer {
    std::vector<float> data;

    OutputBuffer(uint32_t width, uint32_t height);
};

struct AccumulationBuffer {
    std::vector<float> data;
    uint64_t num_samples;

    AccumulationBuffer(uint32_t width, uint32_t height);

    void update_output(OutputBuffer &output);
};