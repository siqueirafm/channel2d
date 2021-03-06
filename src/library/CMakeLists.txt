set (CHANNEL2D_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set (CHANNEL2D_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB HPP_FILES ${CHANNEL2D_INC}/*.hpp)
file (GLOB CPP_FILES ${CHANNEL2D_SRC}/*.cpp)

# Make a static library (omit STATIC if you want user setting to decide)
add_library(Channel2D
  STATIC
    ${HPP_FILES}
    ${CPP_FILES}
  )

# Library need this directories, and its users will need them too.
target_include_directories(Channel2D
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${CHANNEL2D_INC}
  PRIVATE
    ${GLPK_INCLUDE_DIR}
)

# All users of this library will need at list C++ 11
target_compile_features(Channel2D PUBLIC cxx_std_11)

# Set compiler options for MS Visual Studio and everything else
if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(Channel2D PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(Channel2D PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  list (APPEND flags )
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(Channel2D PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(Channel2D PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

source_group(
  TREE "${CHANNEL2D_INC}"
  PREFIX "Header Files"
  FILES ${HPP_FILES})

source_group(
  TREE "${CHANNEL2D_SRC}"
  PREFIX "Source Files"
  FILES ${CPP_FILES})

install(TARGETS Channel2D DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION include)
