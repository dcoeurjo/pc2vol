project(pc2vol)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(WITH_LIBIGL ON)
SET(WITH_CGAL ON)
SET(WITH_GMP ON)
include(dgtal)
include(polyscope)

add_executable(pc2vol pc2vol.cpp)
target_link_libraries(pc2vol polyscope DGtal ${DGTAL_LIBRARIES} igl::core)
