project(libmesh)

file(GLOB SRCS
  "*.h" "*.hh" "*.hpp"
  "*.c" "*.cc" "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} OpenMeshCore)
add_dependencies     (${PROJECT_NAME} OpenMeshCore)

target_link_libraries(${PROJECT_NAME} OpenMeshMore)
add_dependencies     (${PROJECT_NAME} OpenMeshMore)

target_include_directories(${PROJECT_NAME} PUBLIC .)
