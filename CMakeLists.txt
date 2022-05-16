cmake_minimum_required(VERSION 3.0)

project(Nexus)

# add pop count, c++ version and compiler warnings
add_compile_options(-mpopcnt -std=gnu++14 -Wall -pedantic)

add_library(AuxiliaryClasses src/range.cpp src/fmindex.cpp src/nucleotide.cpp src/searchstrategy.cpp src/buildIndexAuxiliary.cpp src/fmpos.cpp src/fmocc.cpp src/search.cpp src/fmindexDBG.cpp src/strainfreemapper.cpp)

#add_compile_definitions(ALPHABET=6)
add_definitions(-DALPHABET=6)

add_executable(columba src/main.cpp   )
target_link_libraries(columba AuxiliaryClasses )

add_executable(nexus src/mainDBG.cpp src/benchmarking.cpp)
target_link_libraries(nexus AuxiliaryClasses )

add_executable(nexusBuild src/buildDBG.cpp)
target_link_libraries(nexusBuild AuxiliaryClasses )

add_executable(createStyles src/createStyles.cpp)

add_executable(visualizeRead src/visualizeRead.cpp)
target_link_libraries(visualizeRead AuxiliaryClasses)

add_executable(visualizePath src/visualizePath.cpp)
target_link_libraries(visualizePath AuxiliaryClasses)


# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# set include path for Google's sparse hash table
find_package(SparseHash)
if (SPARSEHASH_FOUND)
    include_directories(${SPARSEHASH_INCLUDE_DIR})
else (SPARSEHASH_FOUND)
    message(FATAL_ERROR "\nFATAL ERROR: The required Google SparseHash package"
            " could not be found on this system.  Please refer to the Nexus "
            "manual for the Google Sparsehash installation instructions.  If "
            "you installed Google Sparsehash in a non-standard location "
            "(e.g. somewhere in your homedir), you can point cmake to the "
            "installation location as follows: \ncmake "
            "-DSPARSEHASH_INCLUDE_DIR=<path-to-sparsehash>/include .")
endif(SPARSEHASH_FOUND)

install(TARGETS columba DESTINATION bin)
install(TARGETS nexus DESTINATION bin)
install(TARGETS nexusBuild DESTINATION bin)
install(TARGETS createStyles DESTINATION bin)
install(TARGETS visualizePath DESTINATION bin)
install(TARGETS visualizeRead DESTINATION bin)

# add_subdirectory(unittest)

