cmake_minimum_required(VERSION 3.0)

project(CONSTRAINED-DELAUNAY-TRIANGULATION VERSION 1.0 LANGUAGES CXX)

if (WIN32)
    add_definitions(
        -DNOMINMAX
        -D_USE_MATH_DEFINES
        -D_CRT_SECURE_NO_WARNINGS)
endif ()

find_package(OpenMesh REQUIRED)

add_subdirectory(3rdparty/predicates)
add_subdirectory(src/mesh)
add_subdirectory(src/cdt)
add_subdirectory(test/cdt)
