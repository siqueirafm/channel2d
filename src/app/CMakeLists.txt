cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(channel2d-app LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message (STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}")

if (NOT DEFINED CHANNEL2D_INCLUDE_DIR)
  set (CHANNEL2D_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
  message (STATUS "CHANNEL2D_INCLUDE_DIR is ${CHANNEL2D_INCLUDE_DIR}")
endif()

if (NOT DEFINED GLPK_INCLUDE_DIR)
  set (GLPK_INCLUDE_DIR /usr/include)
  message (STATUS "GLPK_INCLUDE_DIR is ${GLPK_INCLUDE_DIR}")
endif()

if (NOT DEFINED CHANNE2D_LIB_DIR)
  set (CHANNEL2D_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
  message (STATUS "CHANNEL2D_LIB_DIR is ${CHANNEL2D_LIB_DIR}")
endif()

if (NOT DEFINED GLPK_LIB_DIR)
  set (GLPK_LIB_DIR /usr/lib)
  message (STATUS "GLPK_LIB_DIR is ${GLPK_LIB_DIR}")
endif()

set (CHANNEL2DAPP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB CPP_FILES ${CHANNEL2DAPP_SRC}/*.cpp)

add_executable(channel2d-app "")

target_compile_features(channel2d-app PRIVATE cxx_std_11)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(channel2d-app PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(channel2d-app PRIVATE /W4 /WX- /MDd)
  endif ()
else ()
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(channel2d-app PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(channel2d-app PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

target_sources(channel2d-app PRIVATE ${CPP_FILES})

target_include_directories(
  channel2d-app
    PRIVATE
       ${CHANNEL2D_INCLUDE_DIR}
       ${GLPK_INCLUDE_DIR}
)

target_link_libraries(
  channel2d-app
    PRIVATE
      -L${CHANNEL2D_LIB_DIR} -lchannel2d
      -L${GLPK_LIB_DIR} -lglpk
)

install(TARGETS channel2d-app DESTINATION bin)
