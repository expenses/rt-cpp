mkdir -p build && \
cd build && \
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
ninja && \
cd ..
