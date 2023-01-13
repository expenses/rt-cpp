mkdir -p build && \
cd build && \
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_COLOR_DIAGNOSTICS=ON && \
ninja && \
cd ..
