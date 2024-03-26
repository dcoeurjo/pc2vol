if (TARGET DGtal)
  return()
endif()

include(FetchContent)

message(STATUS "Fetching DGtal")

SET(BUILD_EXAMPLES OFF)

FetchContent_Declare(
    DGtal
    GIT_REPOSITORY https://github.com/DGtal-team/DGtal.git
    GIT_TAG c3c5b19c561ecfe36756dbfd88dc38a3dd7ecfd2
    )
FetchContent_MakeAvailable(DGtal)

include("${DGtal_BINARY_DIR}/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")
