if (TARGET DGtal)
  return()
endif()

include(FetchContent)

message(STATUS "Fetching DGtal")

SET(BUILD_EXAMPLES OFF)

FetchContent_Declare(
    DGtal
    GIT_REPOSITORY https://github.com/DGtal-team/DGtal.git
    GIT_TAG 640fd9e755acd6228f43708f064e96f5e82ca8dc
    )
FetchContent_MakeAvailable(DGtal)

include("${DGtal_BINARY_DIR}/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")
