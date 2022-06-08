set (CHANNEL2DAPP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB CPP_FILES ${CHANNEL2DAPP_SRC}/*.cpp)

add_executable(Channel2D-App "")

target_compile_features(Channel2D-App PRIVATE cxx_std_11)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(Channel2D-App PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(Channel2D-App PRIVATE /W4 /WX- /MDd)
  endif ()
else ()
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(Channel2D-App PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(Channel2D-App PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

target_sources(Channel2D-App PRIVATE ${CPP_FILES})

target_include_directories(
  Channel2D-App
    PRIVATE
       ${CHANNEL2D_INCLUDE_DIR}
       ${GLPK_INCLUDE_DIR}
)

target_link_libraries(
  Channel2D-App
    PRIVATE
      -L${CHANNEL2D_LIB_DIR} -lChannel2D
      -L${GLPK_LIB_DIR} -lglpk
)

install(TARGETS Channel2D-App DESTINATION bin)
