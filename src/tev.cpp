#include <arpa/inet.h>
#include <netinet/tcp.h>
#include <span>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <sys/socket.h>
#include <unistd.h>

struct TevConnection {
    int socket_fd;

    TevConnection() {
        socket_fd = socket(AF_INET, SOCK_STREAM, 0);

        if (socket_fd < 0) {
            printf("Failed to create socket: %i\n", socket_fd);
            throw std::runtime_error("Socket creation failed");
        }

        struct sockaddr_in ip_address;
        ip_address.sin_family = AF_INET;
        ip_address.sin_port = htons(14158);
        ip_address.sin_addr.s_addr = inet_addr("127.0.0.1");

        auto connection_result = connect(
            socket_fd, (struct sockaddr *)&ip_address, sizeof(ip_address));

        if (connection_result < 0) {
            printf("connecting failed with %i\n", connection_result);
            throw std::runtime_error("Tev connection failed");
        }
    }

    void send_create(const std::string &image_name, uint32_t width,
                     uint32_t height, std::span<std::string> channel_names) {
        uint8_t packet_type = 4;
        uint8_t grab_focus = 0;
        uint32_t num_channels = channel_names.size();
        uint8_t null_byte = 0;

        uint32_t length = 4 + (1 * 2) + image_name.length() + 1 + (4 * 3);
        for (auto &channel_name : channel_names) {
            length += channel_name.length() + 1;
        }

        write(socket_fd, &length, sizeof(length));
        write(socket_fd, &packet_type, sizeof(packet_type));
        write(socket_fd, &grab_focus, sizeof(grab_focus));

        write(socket_fd, image_name.data(), image_name.length());
        write(socket_fd, &null_byte, sizeof(null_byte));

        write(socket_fd, &width, sizeof(width));
        write(socket_fd, &height, sizeof(height));
        write(socket_fd, &num_channels, sizeof(num_channels));

        for (auto &channel_name : channel_names) {
            write(socket_fd, channel_name.data(), channel_name.length());
            write(socket_fd, &null_byte, sizeof(null_byte));
        }
    }

    void send_update(const std::string &image_name, uint32_t x, uint32_t y,
                     uint32_t width, uint32_t height,
                     std::span<std::string> channel_names,
                     std::span<uint64_t> channel_offsets,
                     std::span<uint64_t> channel_strides,
                     std::span<float> data) {
        uint8_t packet_type = 6;
        uint8_t grab_focus = 0;
        uint32_t channel_count = channel_names.size();
        uint8_t null_byte = 0;

        uint32_t length = 4 + (1 * 2) + image_name.length() + 1 + (4 * 5) +
                          channel_offsets.size() * 8 +
                          channel_strides.size() * 8 + data.size() * 4;
        for (auto &channel_name : channel_names) {
            length += channel_name.length() + 1;
        }

        write(socket_fd, &length, sizeof(length));
        write(socket_fd, &packet_type, sizeof(packet_type));
        write(socket_fd, &grab_focus, sizeof(grab_focus));

        write(socket_fd, image_name.data(), image_name.length());
        write(socket_fd, &null_byte, sizeof(null_byte));

        write(socket_fd, &channel_count, sizeof(channel_count));

        for (auto &channel_name : channel_names) {
            write(socket_fd, channel_name.data(), channel_name.length());
            write(socket_fd, &null_byte, sizeof(null_byte));
        }

        write(socket_fd, &x, sizeof(x));
        write(socket_fd, &y, sizeof(y));
        write(socket_fd, &width, sizeof(width));
        write(socket_fd, &height, sizeof(height));

        write(socket_fd, channel_offsets.data(), channel_offsets.size() * 8);
        write(socket_fd, channel_strides.data(), channel_strides.size() * 8);

        write(socket_fd, data.data(), data.size() * 4);
    }
};
