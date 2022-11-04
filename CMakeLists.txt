project(plane-cleaners)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS -O3)

find_package(SDL2 REQUIRED)

add_executable(plane-cleaners src/main.cpp)
target_include_directories(plane-cleaners PUBLIC pcg-cpp/include)
target_link_libraries(plane-cleaners PUBLIC SDL2::SDL2)
