if (TARGET DGtal)
  return()
endif()

include(FetchContent)

message(STATUS "Fetching DGtal")

SET(BUILD_EXAMPLES OFF)

FetchContent_Declare(
    DGtal
    GIT_REPOSITORY https://github.com/DGtal-team/DGtal.git
    GIT_SHALLOW    TRUE
    GIT_TAG bebcc68511d414e6e8a981e666e2b6cc90516135
    )
FetchContent_MakeAvailable(DGtal)

include("${DGtal_BINARY_DIR}/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")
