project(cdt)

file(GLOB SRCS "*.h" "*.hh" "*.hpp" "*.c" "*.cc" "*.cpp")

add_executable(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} libcdt)
add_dependencies     (${PROJECT_NAME} libcdt)