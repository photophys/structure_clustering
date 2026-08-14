cmake_minimum_required(VERSION 3.15...3.29)

project(_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include_directories("boost/boost")

file(GLOB CPP_FILES src/*.cpp)
pybind11_add_module(_core ${CPP_FILES})

install(TARGETS _core DESTINATION structure_clustering)
