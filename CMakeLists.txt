cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(
  Channel2D
  VERSION 1.0
  DESCRIPTION "A library to thread spline curves into polygonal channels"
  LANGUAGES CXX)

message (STATUS "CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")

if (NOT CMAKE_BUILD_TYPE)
  message (STATUS "Setting build type to 'Release' as none was specified")
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message (STATUS "Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")

# enable folders for projects in Visual Studio
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

if (NOT DEFINED CHANNEL2D_INCLUDE_DIR)
  set (CHANNEL2D_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
  message (STATUS "CHANNEL2D_INCLUDE_DIR is ${CHANNEL2D_INCLUDE_DIR}")
endif()

if (NOT DEFINED CHANNE2D_LIB_DIR)
  set (CHANNEL2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
  message (STATUS "CHANNEL2D_LIB_DIR is ${CHANNEL2D_LIB_DIR}")
endif()

if (NOT DEFINED GLPK_INCLUDE_DIR)
  set (GLPK_INCLUDE_DIR /usr/include)
  message (STATUS "GLPK_INCLUDE_DIR is ${GLPK_INCLUDE_DIR}")
endif()

if (NOT DEFINED GLPK_LIB_DIR)
  set (GLPK_LIB_DIR /usr/lib)
  message (STATUS "GLPK_LIB_DIR is ${GLPK_LIB_DIR}")
endif()

if (NOT DEFINED GLPK_LIB_FILENAME)
  set (GLPK_LIB_FILENAME glpk)
  message (STATUS "GLPK_LIB_FILENAME is ${GLPK_LIB_FILENAME}")
endif()

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # ensure -std=c++xx instead of -std=g++xx
  set (CMAKE_CXX_EXTENSIONS OFF)

  # support folders in IDE's
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # indicate documentation build as an option and set it to ON by default
  option(BUILD_DOC "Build documentation" ON)

  # docs only available if this is the main app
  find_package(Doxygen QUIET)
  if (Doxygen_FOUND)
    add_subdirectory(doc)
    install(CODE  "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target doc)")
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/doc)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/latex DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/doc)
  else()
    message (STATUS "Doxygen not found, not building docs")
  endif()
endif()

# The channel2d library code is here
add_subdirectory(src/library)

# The channel2d application code is here
add_subdirectory(src/app)


