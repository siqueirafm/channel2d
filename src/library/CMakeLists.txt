cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(channel2d LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message (STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}")

if (NOT DEFINED GLPK_INCLUDE_DIR)
  set (GLPK_INCLUDE_DIR /usr/include)
  message (STATUS "GLPK_INCLUDE_DIR is ${GLPK_INCLUDE_DIR}")
endif()

if (NOT DEFINED GLPK_LIB_DIR)
  set (GLPK_LIB_DIR /usr/lib)
  message (STATUS "GLPK_LIB_DIR is ${GLPK_LIB_DIR}")
endif()

set (CHANNEL2D_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set (CHANNEL2D_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB HPP_FILES ${CHANNEL2D_INC}/*.hpp)
file (GLOB CPP_FILES ${CHANNEL2D_SRC}/*.cpp)

add_library(channel2d
  STATIC
    ${HPP_FILES}
    ${CPP_FILES}
  )

target_include_directories(channel2d
  PUBLIC 
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${GLPK_INCLUDE_DIR}
    ${CHANNEL2D_INC}
)

target_compile_features(channel2d PRIVATE cxx_std_11)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(channel2d PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(channel2d PRIVATE /W4 /WX- /MDd)
  endif ()
else ()
  list (APPEND flags )
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(channel2d PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(channel2d PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

install(TARGETS channel2d DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION include)
