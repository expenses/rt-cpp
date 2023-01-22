#cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS=-ftime-trace -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_COLOR_DIAGNOSTICS=ON && \
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_COLOR_DIAGNOSTICS=ON && \
cmake --build build
