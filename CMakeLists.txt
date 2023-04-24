cmake_minimum_required (VERSION 2.6)

project (TPFluid)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (common)

# set options
if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# required packages
find_package (glfw3 REQUIRED)

find_package (Eigen REQUIRED) 
include_directories (${EIGEN_INCLUDE_DIR})

find_package (LAPACK REQUIRED)
add_definitions (-DHAVE_LAPACK)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})

find_package (BLAS REQUIRED)
add_definitions (-DHAVE_BLAS)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${BLAS_LIBRARIES})

find_package (OpenMP REQUIRED)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)

find_package (CUDAToolkit REQUIRED)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} CUDA::cudart)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/util)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_compile_definitions(THRUST_IGNORE_CUB_VERSION_CHECK)
add_compile_definitions(_GNU_SOURCE)
add_compile_definitions(STB_IMAGE_WRITE_IMPLEMENTATION)

add_subdirectory (LosTopos)

add_subdirectory (src)

