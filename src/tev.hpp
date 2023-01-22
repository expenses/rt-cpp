#pragma once
#include <cstdint>
#include <span>
#include <string>

struct TevConnection {
    int socket_fd;

    TevConnection();

    void send_create(const std::string &image_name, uint32_t width, uint32_t height,
                     const std::span<std::string> channel_names);

    void send_update(const std::string &image_name, uint32_t x, uint32_t y, uint32_t width, uint32_t height,
                     const std::span<std::string> channel_names, std::span<uint64_t> channel_offsets,
                     std::span<uint64_t> channel_strides, std::span<float> data);
};
