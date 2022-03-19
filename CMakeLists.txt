# ----------------------------------------------------------------------------
# kCluster project
# ----------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.14...3.20)
# Name of project and that it is C/C++.
project (
        refDecontam
        LANGUAGES CXX C
        VERSION 2.0.0
)

set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lstdc++fs -fPIC -lgomp -lrt -fopenmp -pedantic -O3 -Ofast")
# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(lib/kProcessor)

include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/MQF/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/kmerDecoder/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/kmerDecoder/lib/kseq/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/sdsl-lite/include")


set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/kProcessor/include" "${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor")
include_directories(${kProcessor_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/argh")

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/kProcessor/ThirdParty/kmerDecoder/lib/parallel-hashmap")
include_directories("${PHMAP_INCLUDE_DIRS}")

add_library(refDecontam STATIC  ${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor)
set_target_properties(refDecontam PROPERTIES POSITION_INDEPENDENT_CODE 1 CXX_STANDARD 17)
target_link_libraries (refDecontam kProcessor z)
target_include_directories(refDecontam INTERFACE ${PHMAP_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor)

add_executable(decontaminate src/ref_cont.cpp)
target_link_libraries(decontaminate refDecontam kProcessor z)

add_executable(decontam_by_sketches src/decontam_by_sketches.cpp)
target_link_libraries(decontam_by_sketches refDecontam kProcessor z)
