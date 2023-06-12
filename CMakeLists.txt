cmake_minimum_required(VERSION 3.21)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

project(gr_GCMC_clion)

include(CPM.cmake)
include(ExternalProject)

##########

# GOOGLE TEST

##########
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

##########
# JSON
##########

CPMAddPackage(
        NAME nlohmann_json VERSION 3.9.1
        URL https://github.com/nlohmann/json/releases/download/v3.9.1/include.zip
        OPTIONS "JSON_BuildTests OFF"
)

if (nlohmann_json_ADDED)
    add_library(nlohmann_json INTERFACE IMPORTED)
    target_include_directories(nlohmann_json INTERFACE ${nlohmann_json_SOURCE_DIR}/include)
endif()

add_compile_definitions("NLOHMANN_JSON_HPP") # older versions used this macro. Now it's suffixed with "_"


##########
# XRDFILE
##########

ExternalProject_Add(
        project_xdrfile
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
        URL "https://github.com/chemfiles/xdrfile/archive/8935d749e1f43a87221089588d1cc3f37a0354b0.tar.gz"
        URL_HASH SHA256=a5530703fd07a5baadc9ba75d806fe0844d7b3da0e16f5adbb966660a1cd6828
        PATCH_COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/cmake/patches/xdrfile-01.patch
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} xdrfile-static
        UPDATE_DISCONNECTED ON
        CMAKE_ARGS -Wno-dev -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
        LOG_DOWNLOAD ON INSTALL_COMMAND "")

ExternalProject_Get_Property(project_xdrfile source_dir)
ExternalProject_Get_Property(project_xdrfile binary_dir)
set(XdrfileIncludeDir ${source_dir}/include)
add_library(xdrfile STATIC IMPORTED)
set_property(TARGET xdrfile PROPERTY IMPORTED_LOCATION ${binary_dir}/libxdrfile-static.a)
add_dependencies(xdrfile project_xdrfile)
set_target_properties(xdrfile PROPERTIES POSITION_INDEPENDENT_CODE TRUE)


include_directories(${XdrfileIncludeDir})


###################
# PROGRESS TRACKER
###################

ExternalProject_Add(
        project_progresstracker
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
        CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_POSITION_INDEPENDENT_CODE=on
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} progresstracker
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
        URL https://github.com/mlund/progress-cpp/archive/74c33b1eb21417fef9e5fc2b02c7dbe1d533010c.zip
        URL_HASH SHA256=45e2e83a351d44fc1723aecdf1fbf7cee1afc5d44b7190128d8fd6b4437d15b4
)
ExternalProject_Get_Property(project_progresstracker binary_dir)
ExternalProject_Get_Property(project_progresstracker source_dir)
set(ProgressTrackerIncludeDir ${source_dir})
add_library(progresstracker STATIC IMPORTED GLOBAL)
add_dependencies(progresstracker project_progresstracker)
set_property(TARGET progresstracker PROPERTY IMPORTED_LOCATION ${binary_dir}/libprogresstracker.a)



# ========== pre-compiler system headers ==========

option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if (ENABLE_PCH)
    target_precompile_headers(
            gr_GCMC_clion INTERFACE <vector> <string> <map> <set> <algorithm> <utility> <random> <iostream>
            <limits> <memory>
            <nlohmann/json.hpp>)
endif()


# -Rpass-analysis=loop-vectorize -Rpass-missed=loop-vectorize ///I keep this here so that I dont have to google it///

find_package(OpenMP REQUIRED)
find_library(LIB_OPENMP NAMES libiomp5.so)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE   " -O3 -std=c++17 -mavx2 -ffast-math -Rpass-analysis=loop-vectorize -Rpass-missed=loop-vectorize")
set(CMAKE_C_FLAGS_RELEASE   " -O3 -std=c++17 -mavx2 -ffast-math -Rpass-analysis=loop-vectorize -Rpass-missed=loop-vectorize")

#[[set(CMAKE_CXX_FLAGS_DEBUG    " -Og")
set(CMAKE_C_FLAGS_DEBUG   " -Og " )]]

add_executable(gr_GCMC_clion src/main.cpp src/core.h src/TrajectoryAnalysis.cpp src/TrajectoryAnalysis.h src/Readers.cpp
        src/Readers.h src/System.cpp src/System.h src/Energy.cpp src/Energy.h
        src/Grid.cpp src/Grid.h src/VerletList.cpp src/VerletList.h src/vectypes.h src/real.h src/TabelatedLennardJones.h
        src/Coulomb.h src/Systemx.cpp src/Systemx.h src/gmxcomplex.h src/Moves.h src/NeighbourSearcherImplementation.cpp src/NeighbourSearcherImplementation.h)
target_link_libraries(gr_GCMC_clion PRIVATE nlohmann_json xdrfile)


##### THIS IS FOR GOOGLE TEST
enable_testing()

add_executable(
        test_neigbhour_serach
        src/test_neighbour_search.cc
        src/core.h src/TrajectoryAnalysis.cpp src/TrajectoryAnalysis.h src/Readers.cpp
        src/Readers.h src/System.cpp src/System.h src/Energy.cpp src/Energy.h src/NeighbourSearcherImplementation.h src/NeighbourSearcherImplementation.cpp
        src/Grid.cpp src/Grid.h src/VerletList.cpp src/VerletList.h src/vectypes.h src/real.h src/Systemx.cpp src/Systemx.h)
target_link_libraries(
        test_neigbhour_serach
        GTest::gtest_main
        nlohmann_json xdrfile
)

include(GoogleTest)
gtest_discover_tests(test_neigbhour_serach)

#[[
enable_testing()

add_executable(
        test_coulomb
        src/test_coulomb.cc
        src/core.h src/TrajectoryAnalysis.cpp src/TrajectoryAnalysis.h src/Readers.cpp
        src/Readers.h src/System.cpp src/System.h src/Energy.cpp src/Energy.h
        src/Grid.cpp src/Grid.h src/VerletList.cpp src/VerletList.h src/vectypes.h src/real.h src/TabelatedLennardJones.h src/Coulomb.h src/Systemx.cpp src/Systemx.h)
target_link_libraries(
        test_coulomb
        GTest::gtest_main
        nlohmann_json xdrfile
)

include(GoogleTest)
gtest_discover_tests(test_coulomb)]]


