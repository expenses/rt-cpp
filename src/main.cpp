#include <chrono>
#include <iostream>
#include <random>

#include "buffers.cpp"
#include "tev.cpp"

int main() {
  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_real_distribution<float> float_dist(-1.0, 3.0);

  TevConnection connection = TevConnection();

  printf("Started\n");

  std::string image_name = "output";
  std::vector<std::string> channel_names = {"R", "G", "B"};
  std::vector<uint64_t> channel_offsets = {0, 1, 2};
  std::vector<uint64_t> channel_strides = {3, 3, 3};

  uint32_t width = 512;
  uint32_t height = 512;

  auto accum = AccumulationBuffer(width, height);
  auto output = OutputBuffer(width, height);

  connection.send_create(image_name, width, height, channel_names);

  auto start = std::chrono::steady_clock::now();
  bool never_updated = true;

  while (true) {
    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        accum.data[(y * width + x) * 3 + 0] += 0.1 * float_dist(e1);
        accum.data[(y * width + x) * 3 + 1] += 0.2 * float_dist(e1);
        accum.data[(y * width + x) * 3 + 2] += 0.3 * float_dist(e1);
      }
    }

    accum.num_samples += 1;

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (elapsed_seconds.count() > 0.25 || never_updated) {
      printf("accum.num_samples: %lu\n", accum.num_samples);

      start = end;
      never_updated = false;

      accum.update_output(output);

      connection.send_update(image_name, 0, 0, width, height, channel_names,
                             channel_offsets, channel_strides, output.data);
    }
  }
}
